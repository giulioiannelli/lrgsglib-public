from lrgsglib.core import *
from parsers.shared import *
#
description = """
    Computational resources regarding the Signed Laplacian spectrum of 2D 
    lattices.
"""
phelp_L = """
    Size of the square lattice.
"""
phelp_p = """
    Edge flipping probability.
"""
#
DEFAULT_BINSC = 500
DEFAULT_CELL = 'rand'
DEFAULT_EIGMODE = 'scipy'
DEFAULT_GEO = 'squared'
DEFAULT_MODE = 'eigvec_dist'  # Setting default mode to 'eigvec_dist' for contextual relevance
DEFAULT_NAVG = 1000
DEFAULT_PERIOD = 100
DEFAULT_WORKDIR = ''
DEFAULT_HOWMANY = 1
DEFAULT_VERBOSE = False
#
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
    Mode of operation, either 'eigvec_dist' for eigenvector distribution or 'eigval_dist' for eigenvalue distribution | default='{DEFAULT_MODE}'
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
phelp_verbose = f"""
    Increase output verbosity
"""
#
parsDict = {'L': {'help': phelp_L, 'type': int},
            'p': {'help': phelp_p, 'type': float}
            }
#
parsDictOpt = {'mode': {'names': ['-m', '--mode'], 
                        'help': phelp_mode, 
                        'type': str, 
                        'default': DEFAULT_MODE},
                'eigmode': {'names': ['-em', '--eigen_mode'], 
                        'help': phelp_eigMode, 
                        'type': str, 
                        'default': DEFAULT_EIGMODE},
                'geo': {'names': ['-g', '--geo'],
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
                        'default': DEFAULT_WORKDIR}
            }
#
parDA = {'verbose': {'names': ['-v', '--verbose'],
                    'help': phelp_verbose,
                    'action': argparse.BooleanOptionalAction,
                    'default': DEFAULT_VERBOSE}
        }
#
parser = argparse.ArgumentParser(description=description.strip())
# Mandatory arguments
def parsDict_get(var, key):
    return parsDict[var].get(key, None)
for v in parsDict.keys():
    parser.add_argument(v,
        help=parsDict_get(v, 'help'),
        type=parsDict_get(v, 'type'),
    )
# Optional parameters
def parsDictOpt_get(var, key):
    return parsDictOpt[var].get(key, None)
for ov in parsDictOpt.keys():
    parser.add_argument(*parsDictOpt_get(ov, 'names'),
        default=parsDictOpt_get(ov, 'default'),
        help=parsDictOpt_get(ov, 'help'),
        type=parsDictOpt_get(ov, 'type'),
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