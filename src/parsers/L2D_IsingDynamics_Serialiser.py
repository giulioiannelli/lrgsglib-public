from lrgsglib.core import *
#
description = """
    Serialiser for Lattice2D_IsingDynamics.py
"""
L2D_IsingDynamics_progName = "L2D_IsingDynamics"
L2D_IsingDynamics_progNameShrt = "L2DID"
# Default values for the optional parameters
DEFAULT_GEO = 'squared'
DEFAULT_CELL = 'rand'
DEFAULT_INIT_COND = 'ground_state_0'
DEFAULT_RUNLANG = 'C1'
DEFAULT_NAVG = 500
DEFAULT_INSFFX = ''
DEFAULT_OUTSFFX = ''
DEFAULT_NOCLUST = 1
#
DEFAULT_PRINT = False
DEFAULT_EXEC = False
DEFAULT_mMB = 2**10
DEFAULT_MMB = 2**14
# Helpers for argparse descriptions
phelp_print = f"""
    Option to print the output of the Serialiser. 
"""
phelp_exc = f"""
    Option to exec the output of the Serialiser.
"""
phelp_mMB = f"""
    Minimum MB quantity to be allocated for the single process 
"""
phelp_MMB = f"""
    Maximum MB quantity to be allocated for the single process 
"""
phelp_cell = f"""
    Topological defect class: 'rand', 'randXERR', 'randZERR', 'ball_<R>' 
    with type(<R>)=int. 
"""
phelp_geo = f"""
    Geometry of the lattice. 
"""
phelp_navg = f"""
    Number of averages to compute 
"""
phelp_initCond = f"""
    Initial condition for the Ising model
"""
phelp_runlang = f"""
    Language for running the Ising model
"""
phelp_inSuffix = f"""
    Suffix for the input files
"""
phelp_outSuffix = f"""
    Suffix for the output files
"""
phelp_NoClust = f"""
    Number of clusters to compute
"""
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
                    'help': phelp_initCond,
                    'type': str,
                    'default': DEFAULT_INIT_COND},
    'runlang': {'names': ['-rl', '--runlang'],
                'help': phelp_runlang,
                'type': str,
                'default': DEFAULT_RUNLANG},
    'in_suffix': {'names': ['-is', '--in_suffix'],
                    'help': phelp_inSuffix,
                    'type': str,
                    'default': DEFAULT_INSFFX},
    'out_suffix': {'names': ['-os', '--out_suffix'],
                    'help': phelp_outSuffix,
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
            'default': DEFAULT_MMB}
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
parser = argparse.ArgumentParser(description=description)
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