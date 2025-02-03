import numpy as np
# General program arguments
## program helpers
phelp_exc = "Option to exec the output of the Serialiser"
phelp_insuffix = "Suffix for input files"
phelp_mMB = "Minimum MB quantity to be allocated for the single process"
phelp_MMB = "Maximum MB quantity to be allocated for the single process"
phelp_navg = "Number of averages"
phelp_outsuffix = "Suffix for output files"
phelp_print = "Option to print the output of the Serialiser"
phelp_print_chrono = "Print the chronometer"
phelp_remove_files = "Remove the input files after the computation"
phelp_verbose = "Verbose mode"
phelp_workdir = "Working directory"
## default values
DEFAULT_EXEC = False
DEFAULT_INSFFX = ''
DEFAULT_mMB = 2**10
DEFAULT_MMB = 2**14
DEFAULT_NAVG = 500
DEFAULT_OUTSFFX = ''
DEFAULT_PRINT = False
DEFAULT_PRINT_CHRONO = False
DEFAULT_REMOVE_FILES = True
DEFAULT_VERBOSE = False
DEFAULT_WORKDIR = ''
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
# Lattices program arguments
## program helpers
phelp_geo = "Geometry of the lattice"
phelp_side1_list = "List of side1 values for L2D Serializers"
## default values
DEFAULT_GEO = 'squared'
DEFAULT_SIDE1_LIST = [16, 32, 64]
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


