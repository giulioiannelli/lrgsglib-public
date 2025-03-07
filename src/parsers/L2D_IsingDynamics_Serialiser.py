from lrgsglib.core import *
from lrgsglib.config.progargs import *
from parsers.shared import *
#
optional_args_dict = {
    tuple(['-s1', '--side1_list']): {
        'help': phelp_side1_list,
        'type': int,
        'nargs': '+',
        'default': DEFAULT_SIDE1_LIST
    },
    tuple(['-pFT', '--pflip_linsp']): {
        'help': phelp_pflip_linsp,
        'type': parse_multiple_linspace,
        'default': DEFAULT_PFLIP_LINSP
    },
    tuple(['-TT', '--Temp_linsp']): {
        'help': phelp_Temp_linsp,
        'type': parse_multiple_linspace,
        'default': DEFAULT_TEMP_LINSP
    },
    tuple(['-g', '--geometry']): {
        'help': phelp_geo,
        'type': str,
        'default': DEFAULT_GEO
    } ,
    tuple(['-c', '--cell_type']): {
        'help': phelp_cell,
        'type': str,
        'default': DEFAULT_CELL
    },
    tuple(['-n', '--number_of_averages']): {
        'help': phelp_navg,
        'type': int,
        'default': DEFAULT_NAVG
    },
    tuple(['-ic', '--init_cond']): {
        'help': phelp_ic,
        'type': str,
        'default': DEFAULT_INIT_COND
    },
    tuple(['-rl', '--runlang']): {
        'help': phelp_runlang,
        'type': str,
        'default': DEFAULT_RUNLANG
    },
    tuple(['-is', '--in_suffix']): {
        'help': phelp_insuffix,
        'type': str,
        'default': DEFAULT_INSFFX
    },
    tuple(['-os', '--out_suffix']): {
        'help': phelp_outsuffix,
        'type': str,
        'default': DEFAULT_OUTSFFX
    },
    tuple(['-nc', '--NoClust']): {
        'help': phelp_NoClust,
        'type': int,
        'default': DEFAULT_NOCLUST
    },
    tuple(['-mMB', '--slanzarv_minMB']): {
        'help': phelp_mMB,
        'type': int,
        'default': DEFAULT_mMB
    },
    tuple(['-MMB', '--slanzarv_maxMB']): {
        'help': phelp_MMB,
        'type': int,
        'default': DEFAULT_MMB
    },
    tuple(['-ts', '--thermsteps']): {
        'help': phelp_thrmsteps,
        'type': int,
        'default': DEFAULT_THRMSTEPS
    },
    tuple(['-wd', '--workdir']): {
        'help': phelp_workdir,
        'type': str,
        'default': DEFAULT_WORKDIR
    }
}
action_args_dict = {
    tuple(['--exec']): {
        'help': phelp_exc,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_EXEC
    },
    tuple(['--print']): {
        'help': phelp_print,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_PRINT
    }
}
# Setup the argument parser
parser = argparse.ArgumentParser(
    description=L2D_IsingDynamicsSerializer_description, 
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
#
for k,v in {**optional_args_dict, **action_args_dict}.items():
    parser.add_argument(*k, **v)