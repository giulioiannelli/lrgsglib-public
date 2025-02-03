from lrgsglib.core import *
from lrgsglib.config.progargs import *
from parsers.shared import *
#
parserDict_Opt = {
    'side1_list': {'names': ['-s1', '--side1_list'],
                   'help': phelp_side1_list,
                   'type': int,
                   'nargs': '+',
                   'default': DEFAULT_SIDE1_LIST},
    'pflip_linsp': {'names': ['-pFT', '--pflip_linsp'],
                   'help': phelp_pflip_linsp,
                   'type': parse_multiple_linspace,
                   'default': DEFAULT_PFLIP_LINSP},
    'Temp_linsp': {'names': ['-TT', '--Temp_linsp'],
                    'help': phelp_Temp_linsp,
                    'type': parse_multiple_linspace,
                    'default': DEFAULT_TEMP_LINSP},
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
# required arguments


# Optional parameters
for ov in parserDict_Opt.keys():
    parser.add_argument(*parsDict_get(parserDict_Opt, ov, 'names'),
        default=parsDict_get(parserDict_Opt, ov, 'default'),
        help=parsDict_get(parserDict_Opt, ov, 'help')+\
            f"(default: {parsDict_get(parserDict_Opt, ov, 'default')})",
        type=parsDict_get(parserDict_Opt, ov, 'type'),
        nargs=parsDict_get(parserDict_Opt, ov, 'nargs')
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