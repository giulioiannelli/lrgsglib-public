from lrgsglib.core import *
from lrgsglib.config.progargs import *
from parsers.shared import *
#
phelp_L = "Size of the square lattice"
phelp_p = "Edge flipping probability"
phelp_T = "Temperature of the Ising model"
#
parsDict = {
    'L': {'help': phelp_L, 'type': int},
    'p': {'help': phelp_p, 'type': float},
    'T': {'help': phelp_T, 'type': float}
}
#
parsDictOpt = {
    'cell': {
        'names': ['-c', '--cell_type'],
        'help': phelp_cell,
        'type': str,
        'default': DEFAULT_CELL},
    'geo': {
        'names': ['-g', '--geometry'],
        'help': phelp_geo,
        'type': str,
        'default': DEFAULT_GEO},
    'in_suffix': {
        'names': ['-is', '--in_suffix'],
        'help': phelp_insuffix,
        'type': str,
        'default': DEFAULT_INSFFX},
    'init_cond': {
        'names': ['-ic', '--init_cond'],
        'help': phelp_ic,
        'type': str,
        'default': DEFAULT_INIT_COND},
    'navg': {
        'names': ['-n', '--number_of_averages'],
        'help': phelp_navg,
        'type': int,
        'default': DEFAULT_NAVG},
    'NoClust': {
        'names': ['-nc', '--NoClust'],
        'help': phelp_NoClust,
        'type': int,
        'default': DEFAULT_NOCLUST},
    'out_suffix': {
        'names': ['-os', '--out_suffix'],
        'help': phelp_outsuffix,
        'type': str,
        'default': DEFAULT_OUTSFFX},
    'runlang': {
        'names': ['-rl', '--runlang'],
        'help': phelp_runlang,
        'type': str,
        'default': DEFAULT_RUNLANG},
    'thrmsteps': {
        'names': ['-ts', '--thrmsteps'],
        'help': phelp_thrmsteps,
        'type': int,
        'default': DEFAULT_THRMSTEPS},
    'val': {
        'names': ['-V', '--val'],
        'help': phelp_val,
        'type': float,
        'default': DEFAULT_VAL},
    'workdir': {
        'names': ['-wd', '--workdir'],
        'help': phelp_workdir,
        'type': str,
        'default': DEFAULT_WORKDIR},
}
#
parDA = {
    'remove_files': {'names': ['-rf', '--remove_files'],
                     'help': phelp_remove_files,
                     'action': argparse.BooleanOptionalAction,
                     'default': DEFAULT_REMOVE_FILES},
    'verbose': {'names': ['-v', '--verbose'],
                'help': phelp_verbose,
                'action': argparse.BooleanOptionalAction,
                'default': DEFAULT_VERBOSE},
    'print_chrono': {'names': ['-pc', '--print_chrono'],
                     'help': phelp_print_chrono,
                     'action': argparse.BooleanOptionalAction,
                     'default': DEFAULT_PRINT_CHRONO}
}
#
parser = argparse.ArgumentParser(description=L2D_IsingDynamics_description, 
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# Mandatory arguments
for v in parsDict.keys():
    parser.add_argument(v,
        help=parsDict_get(parsDict, v, 'help'),
        type=parsDict_get(parsDict, v, 'type'),
    )
# Optional parameters
for ov in parsDictOpt.keys():
    parser.add_argument(*parsDict_get(parsDictOpt, ov, 'names'),
        default=parsDict_get(parsDictOpt, ov, 'default'),
        help=parsDict_get(parsDictOpt, ov, 'help'),
        type=parsDict_get(parsDictOpt, ov, 'type'),
    )
for ov in parDA.keys():
    parser.add_argument(*parsDict_get(parDA, ov, 'names'),
        default=parsDict_get(parDA, ov, 'default'),
        help=parsDict_get(parDA, ov, 'help'),
        action=parsDict_get(parDA, ov, 'action')
    )