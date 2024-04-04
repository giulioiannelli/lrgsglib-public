from ..shared import *
#
ColorType = Union[Tuple[int, int, int], Tuple[float, float, float], Tuple[int, int, int, int], Tuple[float, float, float, float], str]
# extensions
PDF = ".pdf"
TXT = ".txt"
BIN = ".bin"
MP4 = ".mp4"
PKL = ".pkl"
# paths
PATH_ROOTF = "LRG-Signed"
# default values
LB_PFLIP = .0
UB_PFLIP = 1.

DIR_SRC = "src/"
DIR_PCK = "LRGSG_package/"
DIR_PLT = "plot/"
DIR_DAT = "data/"

DIR_GRAPH = "graphs/"
DIR_ISING = "ising/"
DIR_VOTER = "voter/"
DIR_LRGSG = "lrgsg/"
DIR_PHTRA = "phtra/"
# Signed Graph default values
SG_GRAPH_REPR = 'G'
SG_ERRMSG_PFLIP = f""" pflip must be between {LB_PFLIP} and {UB_PFLIP}, inclusive."""
#
DEFErdosReny_pthabb = "er/"
# 2D Lattice default values
L2D_FBCV = 1.
L2D_PBC = True
L2D_STDFN = ""
L2D_SIDE1 = 32
L2D_SIDE2 = 0
L2D_SGPATH = ""
L2D_ONREP = 'G'
L2D_PHTABB = 'l2d_'
L2D_GEO_TRI = 'triangular'
L2D_GEO_SQR = 'squared'
L2D_GEO_HEX = 'hexagonal'
L2D_GEO_TRI_SHRT = 'tri'
L2D_GEO_SQR_SHRT = 'sqr'
L2D_GEO_HEX_SHRT = 'hex'
L2D_P_C_LIST = [0.146, 0.103, 0.065] 
L2D_WITH_POS = False
# 
L2D_GEO = L2D_GEO_SQR
# 
L2D_GEO_LIST = [L2D_GEO_TRI, 
                        L2D_GEO_SQR, 
                        L2D_GEO_HEX]
L2D_GEO_SHRT_LIST = [L2D_GEO_TRI_SHRT, 
                        L2D_GEO_SQR_SHRT, 
                        L2D_GEO_HEX_SHRT]
L2D_SINGLE_CELL_LIST = ['single', 'singleXERR', 'singleZERR']
L2D_RAND_CELL_LIST = ['rand', 'randXERR', 'randZERR']
#
L2D_P_C_DICT = {g: p for g,p in zip(L2D_GEO_LIST, L2D_P_C_LIST)}
L2D_PATH_DICT = {g: f"{L2D_PHTABB}{g}/" for a,g in 
                        zip(L2D_GEO_SHRT_LIST, L2D_GEO_LIST)}
L2D_GEO_SHRT_DICT = {s: a for a,s in 
                        zip(L2D_GEO_SHRT_LIST, L2D_GEO_LIST)}
L2D_SHRT_GEO_DICT = {a: s for a,s in 
                        zip(L2D_GEO_SHRT_LIST, L2D_GEO_LIST)}

L2D_WARNMSG_GEO = """The selected geometry of the 2D lattice is not available. Setting it to 'squared' for a 2d regular grid."""
L2D_ERRMSG_GEO = """Invalid side value for hexagonal lattice. In order to implement PBC on hexagonal lattice you need to provide an even value for side1 and side2."""

L3D_DIM = (16, 16, 16)
L3D_PBC = True
L3D_FBCV = 1.
L3D_SGPATH = ""
L3D_STDFN = ""
L3D_PHTABB = "l3d_"
L3D_GEO_SC = 'simple_cubic'
L3D_GEO_BCC = 'body_centered'
L3D_GEO_FCC = 'face_centered'
L3D_GEO_SC_SHRT1 = 'cubic'
L3D_GEO_SC_SHRT = 'sc'
L3D_GEO_BCC_SHRT = 'bcc'
L3D_GEO_FCC_SHRT = 'fcc'
L3D_WITH_POS = False
L3D_THETA = np.pi/6
L3D_PHI = np.pi/6
# 
L3D_GEO = L3D_GEO_SC
#
L3D_GEO_LIST = [L3D_GEO_SC, 
                        L3D_GEO_BCC, 
                        L3D_GEO_FCC]
L3D_GEO_SHRT_LIST = [L3D_GEO_SC_SHRT, 
                        L3D_GEO_BCC_SHRT, 
                        L3D_GEO_FCC_SHRT]
L3D_GEO_DICT = {a: g for a,g in 
                 zip(L3D_GEO_SHRT_LIST+L3D_GEO_LIST, L3D_GEO_LIST*2)}
L3D_PATH_DICT = {a: f"{L3D_PHTABB}{g}/" for a,g in 
                 zip(L3D_GEO_SHRT_LIST+L3D_GEO_LIST, L3D_GEO_LIST*2)}

DEFAULT_RECURSION_LIMIT = 1024**2
DEFAULT_ENTROPY_STEPS = 1000
DEFAULT_ENTROPY_LEXPONENT = -3
DEFAULT_ENTROPY_HEXPONENT = 6
DEFAULT_ISING_NOCLUST = 5
DEFAULT_MAX_DIGITS_ROUND_SIGFIG = 18
DEFAULT_NUNMBER_AVERAGES = 100
DEFAULT_SPIKE_THRESHOLD = 0.05
DEFAULT_MAX_THRESHOLD = 2 * DEFAULT_SPIKE_THRESHOLD

# default lists
DEFAULT_ERGRAPH_ABBRV = "randGraph"