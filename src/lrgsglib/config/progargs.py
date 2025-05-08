from parsers.shared import *
import numpy as np
# General program arguments
## program helpers
phelp_binsc = "Number of bins for the distribution sampling"
phelp_data_save_freq = "Frequency of the saved data"
phelp_eigMode = "Library for spectral computations"
phelp_exc = "Option to execute the output of the Serialiser"
phelp_freq = "Frequency of the saved spin output"
phelp_insuffix = "Suffix for input files"
phelp_L = "Main side of the lattice"
phelp_mMB = "Minimum MB quantity to be allocated for the single process"
phelp_MMB = "Maximum MB quantity to be allocated for the single process"
phelp_moretime = "Time allocation for the slanzarv job"
phelp_navg = "Number of averages"
phelp_nomail = "Disable slanzarv email notifications"
phelp_outsuffix = "Suffix for output files"
phelp_p = "Edge flipping probability"
phelp_print = "Option to print the output of the Serialiser"
phelp_print_chrono = "Print the chronometer"
phelp_remove_files = "Remove the input files after the computation"
phelp_T = "Temperature of the Ising model"
phelp_short = "Run slanzarv job in short mode"
phelp_slanzarv_id = "ID for the slanzarv job"
phelp_verbose = "Verbose mode"
phelp_val = "Value for the clusters"
phelp_workdir = "Working directory"
## default values
DEFAULT_BINSC = 500
DEFAULT_EIGMODE = 'scipy'
DEFAULT_EXEC = False
DEFAULT_FREQ = 2
DEFAULT_INSFFX = ''
DEFAULT_mMB = 2**11
DEFAULT_MMB = 2**11
DEFAULT_MORETIME = 0
DEFAULT_NAVG = 500
DEFAULT_NOMAIL = True
DEFAULT_OUTSFFX = ''
DEFAULT_PRINT = False
DEFAULT_PRINT_CHRONO = False
DEFAULT_REMOVE_FILES = True
DEFAULT_SHORT = False
DEFAULT_SLANZARV_ID = ""
DEFAULT_VAL = "=1"
DEFAULT_VERBOSE = False
DEFAULT_WORKDIR = ''
## default values depending on default arguments
DEFAULT_DATA_SAVE_FREQ = DEFAULT_NAVG // 10
# SignedGraph program arguments
## program helpers
phelp_cell = "Topological defect class: 'rand', 'randXERR', 'randZERR', \
    'ball_<R>' with type(<R>)=int"
phelp_NoClust = "Number of clusters to compute"
phelp_pflip_linsp = "Tuple for linspace (e.g., '(0.15, 0.35, 3)')"
phelp_Temp_linsp = "Tuple for linspace (e.g., '(0.1, 1, 10);(1, 2.5, 20);(2.5, 5, 3)')"
## default values
DEFAULT_CELL = 'rand'
DEFAULT_NOCLUST = 1
DEFAULT_PFLIP_LINSP = np.linspace(0.01, 0.5, 10)
DEFAULT_TEMP_LINSP = np.linspace(0.1, 2.2, 10)
# Program arguments
action_args_dict = {
    tuple(['-rf', '--remove_files']): {
        'help': phelp_remove_files,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_REMOVE_FILES},
    tuple(['-v', '--verbose']): {
        'help': phelp_verbose,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_VERBOSE},
    tuple(['-pc', '--print_chrono']): {
        'help': phelp_print_chrono,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_PRINT_CHRONO}
}
# Serializers program arguments
Serializer_action_args_dict = {
    tuple(['--exec']): {
        'help': phelp_exc,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_EXEC
    },
    tuple(['--print']): {
        'help': phelp_print,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_PRINT
    },
    tuple(['--nomail']): {
        'help': phelp_nomail,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_NOMAIL
    },
    tuple(['--short']): {
        'help': phelp_short,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_SHORT
    },
}
# Lattices program arguments
## program helpers
phelp_geo = "Geometry of the lattice"
phelp_side1_list = "List of side1 values for L2D Serializers"
## default values
DEFAULT_GEO = 'squared'
DEFAULT_SIDE1_LIST = [16, 32, 64]
## args parsers dict
L2D_args = {
    'L': {
        'help': phelp_L, 
        'type': int
    },
    'p': {
        'help': phelp_p,
        'type': float
    }
}
L2D_opt_args = {
    tuple(['-c', '--cell_type']): {
        'help': phelp_cell,
        'type': str,
        'default': DEFAULT_CELL
    },
    tuple(['-g', '--geometry']): {
        'help': phelp_geo,
        'type': str,
        'default': DEFAULT_GEO
    },
    tuple(['-na', '--number_of_averages']): {
        'help': phelp_navg,
        'type': int,
        'default': DEFAULT_NAVG
    },
    tuple(['-wd', '--workdir']): {
        'help': phelp_workdir,
        'type': str,
        'default': DEFAULT_WORKDIR
    }
}
# Signed Laplacian Spectra program arguments
## program helpers
phelp_howmany_eigs = "Number of eigenvalues to compute"
phelp_l2dsspect_mode = "Mode of operation, either 'eigvec_dist' for \
    eigenvector distribution or 'eigval_dist' for eigenvalue distribution \
    or 'eigvals' for eigenvalues"
## default values
DEFAULT_HOWMANY_EIGS = 1
DEFAULT_L2DSSPECT_MODE = 'eigvec_dist'
## names and descriptions
L2D_SlaplSpect_progName = 'L2D_SlaplSpect'
L2D_SlaplSpect_progNameShrt = 'L2DSS'
L2D_SlaplSpect_description = f"""
    Computational resourses regarding the Signed Laplacian spectrum of 2D 
    lattices: {L2D_SlaplSpect_progName}.py
"""
## arg parsers dict
L2D_SlaplSpect_args = {**L2D_args}
L2D_SlaplSpect_optional_args_dict = {
    tuple(['-bc', '--bins_count']): {
        'help': phelp_binsc,
        'type': int,
        'default': DEFAULT_BINSC
    },
    tuple(['-em', '--eigen_mode']): {
        'help': phelp_eigMode,
        'type': str,
        'default': DEFAULT_EIGMODE
    },
    tuple(['-hm', '--howmany']): {
        'help': phelp_howmany_eigs,
        'type': int,
        'default': DEFAULT_HOWMANY_EIGS
    },
    tuple(['-m', '--mode']): {
        'help': phelp_l2dsspect_mode,
        'type': str,
        'default': DEFAULT_L2DSSPECT_MODE
    },
    tuple(['-prd', '--period']): {
        'help': phelp_data_save_freq,
        'type': int,
        'default': DEFAULT_DATA_SAVE_FREQ
    }
}
L2D_SlaplSpect_action_args_dict = {}
# IsingDynamics program arguments
## program helpers
phelp_ic = "Initial condition for the Ising model"
phelp_runlang = "Language for running the Ising model"
phelp_thrmsteps = "Number of thermalization steps"
## default values
DEFAULT_INIT_COND = 'ground_state_0'
DEFAULT_RUNLANG = 'C1'
DEFAULT_THRMSTEPS = 20
## names and descriptions
L2D_IsingDynamics_progName = 'L2D_IsingDynamics'
L2D_IsingDynamics_progNameShrt = 'L2DID'
L2D_IsingDynamics_description = f"""
    Computational resourses regarding the Ising Dynamics of 2D 
    lattices: {L2D_IsingDynamics_progName}.py
"""
L2D_IsingDynamicsSerializer_description = f"""
    Serialiser for {L2D_IsingDynamics_progName}.py
"""
## arg parsers dict
L2D_IsingDynamics_args = {
    **L2D_args, 
    'T': {
        'help': phelp_T,
        'type': float
    }
}
L2D_IsingDynamics_action_args_dict = {**action_args_dict}
L2D_IsingDynamics_optional_args_dict = {
    tuple(['-fq', '--freq']): {
        'help': phelp_freq,
        'type': int,
        'default': DEFAULT_FREQ
    },
    tuple(['-ic', '--init_cond']): {
        'help': phelp_ic,
        'type': str,
        'default': DEFAULT_INIT_COND
    },
    tuple(['-is', '--in_suffix']): {
        'help': phelp_insuffix,
        'type': str,
        'default': DEFAULT_INSFFX
    },
    tuple(['-nc', '--NoClust']): {
        'help': phelp_NoClust,
        'type': int,
        'default': DEFAULT_NOCLUST
    },
    tuple(['-os', '--out_suffix']): {
        'help': phelp_outsuffix,
        'type': str,
        'default': DEFAULT_OUTSFFX
    },
    tuple(['-rl', '--runlang']): {
        'help': phelp_runlang,
        'type': str,
        'default': DEFAULT_RUNLANG
    },
    tuple(['-ts', '--thrmsteps']): {
        'help': phelp_thrmsteps,
        'type': int,
        'default': DEFAULT_THRMSTEPS
    },
    tuple(['-vl', '--val']): {
        'help': phelp_val,
        'type': str,
        'default': DEFAULT_VAL
    }
}
L2D_IsingDynamicsSerializer_optional_args_dict = {
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
    tuple(['--moretime']): {
        'help': phelp_moretime,
        'type': int,
        'default': DEFAULT_MORETIME
    },
    tuple(['--slanzarv_id']): {
        'help': phelp_slanzarv_id,
        'type': str,
        'default': DEFAULT_SLANZARV_ID
    },
}
L2D_IsingDynamicsSerializer_action_args_dict = {**Serializer_action_args_dict}
