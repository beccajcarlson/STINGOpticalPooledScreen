import sys
from functools import partial
import ops.annotate
import ops.firesnake
from ops.firesnake import Snake
import ops.io
import pandas as pd


WELLS, TILES = ops.firesnake.load_well_tile_list(config['WELL_TILE_LIST'],include=config['INCLUDE_WELLS_TILES'])
TILES = [str(item).zfill(4) for item in TILES]

# display options for saved .tif files (view in ImageJ)
channels = ('DAPI', 'CELLS')
LUTS = [getattr(ops.io, config['LUTS'][x]) for x in channels]

# set paramspaces if a paramsearch mode is selected
if config['MODE'] == 'paramsearch_segmentation':
    (config,
        nuclei_segmentation_paramspace,
        cell_segmentation_paramspace) = ops.firesnake.initialize_paramsearch(config)
elif config['MODE'] == 'paramsearch_read-calling':
    config,read_calling_paramspace = ops.firesnake.initialize_paramsearch(config)
elif config['MODE']!='process':
    raise ValueError(f'MODE="{config["MODE"]}" not recognized, use either "process" or "paramsearch"')
else:
    if any(map(lambda x: isinstance(x,list),[config['THRESHOLD_DAPI'],config['THRESHOLD_CELL']])):
        raise ValueError('Thresholds cannot be lists for MODE="process"')
    if isinstance(config['NUCLEUS_AREA'][0],list):
        raise ValueError('NUCLEUS_AREA cannot be a list of lists for MODE="process"')

# naming convention for input and processed files
input_files_pheno = partial(ops.firesnake.input_files_nocycles_alt,
                      directory=config['INPUT_DIRECTORY'])

input_files = partial(ops.firesnake.input_files,
                      magnification=config['MAGNIFICATION'],
                      directory=config['INPUT_DIRECTORY'])

processed_input = partial(ops.firesnake.processed_file,
                         directory=config['PROCESS_DIRECTORY'],
                         )

processed_output = partial(ops.firesnake.processed_file,
                         directory=config['PROCESS_DIRECTORY'],
                         temp_tags=config['TEMP_TAGS'],
                         )

rule all:
    input:
        # request individual files or list of files
        [expand(processed_input(x), zip, well=WELLS, tile=TILES)
            for x in config['REQUESTED_TAGS']],
        [config['PROCESS_DIRECTORY'] + '/' + x for x in config['REQUESTED_FILES']],

rule illum_corr:
    input:
        input_files_pheno(config['PHENOTYPE_CYCLE'][0])
    output:
        processed_output('phenotype_corr.tif')
    run:
        Snake.find_corr(output=output, data=input[0],smooth=config['CORR_FILE'],
            autoscale=config['AUTOSCALE_PHENOTYPE'])


rule segment:
    input:
        processed_output('phenotype_corr.tif'),
    output:
        processed_output('nuclei.tif'),
        processed_output('cells.tif'),
    run:
        if config['SEGMENT_METHOD'] == 'cell_2019':
            Snake.segment_cell_2019(
                output=output, 
                data=input[0],
                nuclei_threshold=config['THRESHOLD_DAPI'],
                nuclei_area_min=config['NUCLEUS_AREA'][0],
                nuclei_area_max=config['NUCLEUS_AREA'][1],
                cell_threshold=config['THRESHOLD_CELL'],
            )
        elif config['SEGMENT_METHOD'] == 'cell_2019_select_channels':
            Snake.segment_cell_2019_select_channels(
                output=output,
                data=input[0],
                nuclei_threshold=config['THRESHOLD_DAPI'],
                nuclei_area_min=config['NUCLEUS_AREA'][0],
                nuclei_area_max=config['NUCLEUS_AREA'][1],
                cell_threshold=config['THRESHOLD_CELL'],
                chstart=1, chend=2)
        elif config['SEGMENT_METHOD'] == 'cellpose':
            # last cycle
            cycle = config['CELLPOSE']['CYTO_CYCLE']
            data = ops.io.read_stack(input[0])[cycle]
            Snake.segment_cellpose(
                output=output, 
                data=data, 
                dapi_index=0, 
                cyto_index=config['CELLPOSE']['CYTO_CHANNEL'],
                diameter=config['CELLPOSE']['DIAMETER'],
                )
        else:
            error = ('config entry SEGMENT_METHOD must be "cell_2019" or "cellpose", '
                     f'not {config["SEGMENT_METHOD"]}')
            raise ValueError(error)

rule find_cytoplasm:
    input:
        processed_input('nuclei.tif'),
        processed_input('cells.tif'),
    output:
        processed_output('cytoplasm.tif')
    run:
        Snake.find_cytoplasm(output=output, nuclei=input[0],
            cells=input[1])


rule pheno_ch0:
    input:
        processed_input('phenotype_corr.tif'),
        processed_input('nuclei.tif'),
        processed_input('cells.tif'),
        processed_input('cytoplasm.tif'),
    output:
        processed_output('dapi.csv')
    run:
        Snake.extract_phenotype_extended_channel(output=output, data_phenotype=input[0], nuclei=input[1], cells=input[2], cytoplasm=input[3], 
        	channel = 0, corrchannel1 = 1, wildcards=wildcards, chnames = ['DAPI','PSTING'])


rule pheno_ch1:
    input:
        processed_input('phenotype_corr.tif'),
        processed_input('nuclei.tif'),
        processed_input('cells.tif'),
        processed_input('cytoplasm.tif'),
    output:
        processed_output('psting.csv')
    run:
        Snake.extract_phenotype_extended_channel(output=output, data_phenotype=input[0], nuclei=input[1], cells=input[2], cytoplasm=input[3],
                channel = 1, wildcards=wildcards, chnames = ['DAPI','PSTING'])


rule pheno_morphology:
    input:
        processed_input('phenotype_corr.tif'),
        processed_input('nuclei.tif'),
        processed_input('cells.tif'),
        processed_input('cytoplasm.tif'),
    output:
        processed_output('morph.csv')
    run:
        Snake.extract_phenotype_extended_morphology(output=output, data_phenotype=input[0], nuclei=input[1], cells=input[2], cytoplasm=input[3], wildcards=wildcards)

rule call_min_pheno:
     input:
        processed_input('cells.tif'),
        processed_input('nuclei.tif'),
     output:
        processed_output('coords.csv')
     run:
        Snake.extract_phenotype_minimal(output=output, data_phenotype=input[0], nuclei=input[1], wildcards=wildcards)

if config['MODE'] == 'paramsearch_segmentation':
    rule segment_nuclei_paramsearch:
        input:
            processed_output('phenotype_corr.tif'),
            #input_files(config['SBS_INPUT_TAG'], SBS_CYCLES[0]),
            #not used, just here to change the order of rule execution
            #processed_output('phenotype_aligned.tif'),
        output:
            processed_output(f'nuclei.{nuclei_segmentation_paramspace.wildcard_pattern}.tif')
        params:
            nuclei_segmentation = nuclei_segmentation_paramspace.instance
        run:
            Snake.segment_nuclei(output=output, data=input[0],
                threshold=params.nuclei_segmentation['THRESHOLD_DAPI'],
                area_min=params.nuclei_segmentation['NUCLEUS_AREA_MIN'],
                area_max=params.nuclei_segmentation['NUCLEUS_AREA_MAX'])

    rule segment_cells_paramsearch:
        input:
            processed_output('phenotype_corr.tif'),
            #input_files(config['SBS_INPUT_TAG'], SBS_CYCLES[0]),
            processed_input(f'nuclei.{nuclei_segmentation_paramspace.wildcard_pattern}.tif')
        output:
            processed_output(f'cells.{nuclei_segmentation_paramspace.wildcard_pattern}.{cell_segmentation_paramspace.wildcard_pattern}.tif')
        params:
            cell_segmentation = cell_segmentation_paramspace.instance
        run:
            Snake.segment_cells_select_channels(output=output,
                data=input[0], nuclei=input[1], chstart=1, chend=2, cell_threshold=params.cell_segmentation['THRESHOLD_CELL'])

    rule segment_paramsearch_summary:
        input:
            data = processed_output('phenotype_corr.tif'),
            #data = input_files(config['SBS_INPUT_TAG'], SBS_CYCLES[0]),
            segmentations = [processed_input(f'nuclei.{nuclei_segmentation_paramspace.wildcard_pattern}.tif')]+
            [processed_input(f'cells.{nuclei_segmentation_paramspace.wildcard_pattern}.'
                f'{cell_segmentation_instance}.tif')
                for cell_segmentation_instance in cell_segmentation_paramspace.instance_patterns
                ]
        output:
            processed_output(f'segmentation_summary.{nuclei_segmentation_paramspace.wildcard_pattern}.'
                f'{"_".join(cell_segmentation_paramspace.instance_patterns)}.tif')
        run:
            Snake.summarize_paramsearch_segmentation(output=output,data=input[0],segmentations=input.segmentations,
                luts=LUTS[:2]+[ops.io.GLASBEY,]*len(input.segmentations), autoscale=config['AUTOSCALE_PHENOTYPE']
                )
