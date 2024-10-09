from lrgsglib.core import *
#
ER_IsingDynamics_progName = "ER_IsingDynamics"
ER_IsingDynamics_progNameShrt = "ERID"
#
description = f"""
    Serialiser for {ER_IsingDynamics_progName}.py
        - Computes the dynamics of the Ising model on an Erdos-Renyi graph with
            a given average degree and flipping different configuration of 
            edges.
        - CODE ID: {ER_IsingDynamics_progNameShrt}
"""
# Default values for the optional parameters
DEFAULT_CELL = 'rand'
DEFAULT_INIT_COND = 'ground_state_0'
DEFAULT_INSFFX = ''
DEFAULT_OUTSFFX = ''
DEFAULT_RUNLANG = 'C1'
#
DEFAULT_K = 10
DEFAULT_NAVG = 500
DEFAULT_NOCLUST = 1
#
DEFAULT_PRINT = False
DEFAULT_EXEC = False
#
DEFAULT_mMB = 2**10
DEFAULT_MMB = 2**14
#
# Helpers for argparse descriptions
phelp_K = f"""
    Average degree of the Erdos-Renyi graph. | default={DEFAULT_K}
"""
phelp_print = f"""
    Option to print the output of the Serialiser. | default={DEFAULT_PRINT}
"""
phelp_exc = f"""
    Option to exec the output of the Serialiser. | default={DEFAULT_EXEC}
"""
phelp_mMB = f"""
    Minimum MB quantity to be allocated for the single process | 
    default={DEFAULT_mMB}
"""
phelp_MMB = f"""
    Maximum MB quantity to be allocated for the single process | 
    default={DEFAULT_MMB}
"""
phelp_cell = f"""
    Topological defect class: 'rand', 'randXERR', 'randZERR'. | 
    default='{DEFAULT_CELL}'
"""
phelp_navg = f"""
    Number of averages to compute | default={DEFAULT_NAVG}
"""
phelp_initCond = f"""
    Initial condition for the Ising model | default='{DEFAULT_INIT_COND}'
"""
phelp_runlang = f"""
    Language for running the Ising model | default='{DEFAULT_RUNLANG}'
"""
phelp_inSuffix = f"""
    Suffix for the input files | default='{DEFAULT_INSFFX}'
"""
phelp_outSuffix = f"""
    Suffix for the output files | default='{DEFAULT_OUTSFFX}'
"""
phelp_NoClust = f"""
    Number of clusters to compute | default={DEFAULT_NOCLUST} 
"""
#
parDO = {
    'K': {'names': ['-k', '--average_degree'],
            'help': phelp_K,
            'type': int,
            'default': DEFAULT_K},
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