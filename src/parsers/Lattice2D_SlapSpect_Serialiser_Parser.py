from LRGSG_package.LRGSG import *
#
description = """
    Serialiser for Lattice2D_SlaplSpect.py
"""
# Default values for the optional parameters
DEFAULT_BINSC = 500
DEFAULT_EIGMODE = 'scipy'
DEFAULT_NAVG = 10**4
DEFAULT_PERIOD = DEFAULT_NAVG//20
DEFAULT_WORKDIR = ''
DEFAULT_HOWMANY = 1
DEFAULT_GEO = 'squared'
DEFAULT_CELL = 'rand'
DEFAULT_MODE = 'slanzarv_eigDistr'
#
DEFAULT_PRINT = False
DEFAULT_EXEC = False
DEFAULT_mMB = 2**10
DEFAULT_MMB = 2**14
# Helpers for argparse descriptions
phelp_print = f"""
    Option to print the output of the Serialiser. | default={DEFAULT_PRINT}
"""
phelp_exc = f"""
    Option to exec the output of the Serialiser. | default={DEFAULT_EXEC}
"""
phelp_mMB = f"""
    Minimum MB quantity to be allocated for the single process 
    | default={DEFAULT_mMB}
"""
phelp_MMB = f"""
    Maximum MB quantity to be allocated for the single process 
    | default={DEFAULT_MMB}
"""
phelp_binsc = f"""
    Number of bins for the distribution sampling | default={DEFAULT_BINSC}
"""
phelp_cell = f"""
    Topological defect class: 'rand', 'randXERR', 
    'randZERR', 'ball_<R>' with type(<R>)=int. | default='{DEFAULT_CELL}'
"""
phelp_eigMode = f"""
    Spectral computations (numpy/scipy) | default='{DEFAULT_EIGMODE}'
"""
phelp_geo = f"""
    Geometry of the lattice. | default='{DEFAULT_GEO}'
"""
phelp_howmany = f"""
    Number of eigenvalues to compute | default={DEFAULT_HOWMANY}
"""
phelp_mode = f"""
    Mode, (eigDistr, slanzarv_eigDistr) | default='{DEFAULT_MODE}'
"""
phelp_navg = f"""
    Number of averages to compute | default={DEFAULT_NAVG}
"""
phelp_period = f"""
    Period for saving the data | default={DEFAULT_PERIOD}
"""
phelp_workDir = f"""
    Working directory | default='{DEFAULT_WORKDIR}'
"""
#
parsDictOpt = {'mode': {'names': ['-m', '--mode'], 
                        'help': phelp_mode, 
                        'type': str, 
                        'default': DEFAULT_MODE},
                'eigmode': {'names': ['-em', '--eigen_mode'], 
                        'help': phelp_eigMode, 
                        'type': str, 
                        'default': DEFAULT_EIGMODE},
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
                'period': {'names': ['-prd', '--period'],
                        'help': phelp_period,
                        'type': int,
                        'default': DEFAULT_PERIOD},
                'bins_count': {'names': ['-bc', '--bins_count'],
                        'help': phelp_binsc,
                        'type': int,
                        'default': DEFAULT_BINSC},
                'howmany': {'names': ['-hm', '--howmany'],
                        'help': phelp_howmany,
                        'type': int,
                        'default': DEFAULT_HOWMANY},
                'workDir': {'names': ['-wd', '--workDir'],
                        'help': phelp_workDir,
                        'type': str,
                        'default': DEFAULT_WORKDIR},
                'mMB': {'names': ['-mMB', '--slanzarv_minMB'],               
                        'help': phelp_mMB,
                        'type': int,
                        'default': DEFAULT_mMB},
                'MMB': {'names': ['-MMB', '--slanzarv_maxMB'],
                        'help': phelp_MMB,
                        'type': int,
                        'default': DEFAULT_MMB}  
                        }
parsDictAct = {'exec': {'names': ['-e', '--exec'],
                        'help': phelp_exc,
                        'action': argparse.BooleanOptionalAction,
                        'default': DEFAULT_EXEC},
                'print': {'names': ['-p', '--print'],                 
                        'help': phelp_print,
                        'action': argparse.BooleanOptionalAction,
                        'default': DEFAULT_PRINT}
                        }
Lattice2D_SlaplSpect_progName = "Lattice2D_SlaplSpect"
Lattice2D_SlaplSpect_progNameShrt = "L2DSS"
# Setup the argument parser
parser = argparse.ArgumentParser(description=description)
# Optional parameters
def parsDict_get(var, key):
    return parsDictOpt[var].get(key, None)
for ov in parsDictOpt.keys():
    parser.add_argument(*parsDict_get(ov, 'names'),
        default=parsDict_get(ov, 'default'),
        help=parsDict_get(ov, 'help'),
        type=parsDict_get(ov, 'type'),
    )
# Optional parameters
def parsDict_get(var, key):
    return parsDictAct[var].get(key, None)
for ov in parsDictAct.keys():
    parser.add_argument(*parsDict_get(ov, 'names'),
        default=parsDict_get(ov, 'default'),
        help=parsDict_get(ov, 'help'),
        action=parsDict_get(ov, 'action')
    )