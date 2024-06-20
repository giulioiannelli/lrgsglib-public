from lrgsglib.core import *
#
description = """
    Computational resourses regarding the Ising Dynamics of 2D 
    lattices.
"""
phelp_L = """
    Size of the square lattice.
"""
phelp_p = """
    Edge flipping probability.
"""
phelp_T = """
    Temperature of the Ising model.
"""
#
DEFAULT_GEO = 'squared'
DEFAULT_CELL = 'rand'
DEFAULT_INIT_COND = 'ground_state_0'
DEFAULT_RUNLANG = 'C1'
DEFAULT_NAVG = 500
DEFAULT_INSFFX = ''
DEFAULT_OUTSFFX = ''
DEFAULT_NOCLUST = 1
#
phelp_cell = f"""
    Topological defect class: 'rand', 'randXERR', 
    'randZERR', 'ball_<R>' with type(<R>)=int. | default='{DEFAULT_CELL}'
"""
phelp_geo = f"""
    Geometry of the lattice. | default='{DEFAULT_GEO}'
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
parsDict = {
    'L': {'help': phelp_L, 'type': int},
    'p': {'help': phelp_p, 'type': float},
    'T': {'help': phelp_T, 'type': float}
}
#
parsDictOpt = {
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
                'default': DEFAULT_NOCLUST}
}
#

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
def parsDict_get(var, key):
    return parsDictOpt[var].get(key, None)
for ov in parsDictOpt.keys():
    parser.add_argument(*parsDict_get(ov, 'names'),
        default=parsDict_get(ov, 'default'),
        help=parsDict_get(ov, 'help'),
        type=parsDict_get(ov, 'type'),
    )