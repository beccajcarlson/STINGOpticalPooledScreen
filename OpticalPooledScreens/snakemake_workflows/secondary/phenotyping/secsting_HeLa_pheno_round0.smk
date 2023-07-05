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
            #display_ranges=DISPLAY_RANGES, luts=LUTS)
