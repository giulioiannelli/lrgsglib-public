from lrgsglib.core import *
from lrgsglib.config.progargs import *
#
parDO = {
    'geo': {'names': ['-g', '--geometry'],
            'help': phelp_geo,
            'type': str,
            'default': DEFAULT_GEO} ,
    'cell': {'names': ['-c', '--cell_type'],
            'help': phelp_cell,
            'type': str,
            'default': DEFAULT_CELL},
    'navg': {'names': ['-n', '--number_of_averages'],
            'help': phelp_navg,
            'type': int,
            'default': DEFAULT_NAVG},
    'init_cond': {'names': ['-ic', '--init_cond'],
                    'help': phelp_ic,
                    'type': str,
                    'default': DEFAULT_INIT_COND},
    'runlang': {'names': ['-rl', '--runlang'],
                'help': phelp_runlang,
                'type': str,
                'default': DEFAULT_RUNLANG},
    'in_suffix': {'names': ['-is', '--in_suffix'],
                    'help': phelp_insuffix,
                    'type': str,
                    'default': DEFAULT_INSFFX},
    'out_suffix': {'names': ['-os', '--out_suffix'],
                    'help': phelp_outsuffix,
                    'type': str,
                    'default': DEFAULT_OUTSFFX},
    'NoClust': {'names': ['-nc', '--NoClust'],
                'help': phelp_NoClust,
                'type': int,
                'default': DEFAULT_NOCLUST},
    'mMB': {'names': ['-mMB', '--slanzarv_minMB'],               
            'help': phelp_mMB,
            'type': int,
            'default': DEFAULT_mMB},
    'MMB': {'names': ['-MMB', '--slanzarv_maxMB'],
            'help': phelp_MMB,
            'type': int,
            'default': DEFAULT_MMB},
    'thermsteps': {'names': ['-ts', '--thermsteps'],
                    'help': phelp_thrmsteps,
                    'type': int,
                    'default': DEFAULT_THRMSTEPS},
    'workdir': {'names': ['-wd', '--workdir'],
                'help': phelp_workdir,
                'type': str,
                'default': DEFAULT_WORKDIR}
}
parDA = {'exec': {'names': ['-e', '--exec'],
                        'help': phelp_exc,
                        'action': argparse.BooleanOptionalAction,
                        'default': DEFAULT_EXEC},
                'print': {'names': ['-p', '--print'],                 
                        'help': phelp_print,
                        'action': argparse.BooleanOptionalAction,
                        'default': DEFAULT_PRINT}
                        }

# Setup the argument parser
parser = argparse.ArgumentParser(description=L2D_IsingDynamicsSerializer_description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# Optional parameters
def parsDict_get(var, key):
    return parDO[var].get(key, None)
for ov in parDO.keys():
    parser.add_argument(*parsDict_get(ov, 'names'),
        default=parsDict_get(ov, 'default'),
        help=parsDict_get(ov, 'help')+\
            f"(default: {parsDict_get(ov, 'default')})",
        type=parsDict_get(ov, 'type'),
    )
# Optional parameters
def parsDict_get(var, key):
    return parDA[var].get(key, None)
for ov in parDA.keys():
    parser.add_argument(*parsDict_get(ov, 'names'),
        default=parsDict_get(ov, 'default'),
        help=parsDict_get(ov, 'help'),
        action=parsDict_get(ov, 'action')
    )