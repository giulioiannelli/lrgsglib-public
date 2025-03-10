from ..shared import *
# types
ColorType = Union[
    Tuple[int, int, int],
    Tuple[float, float, float],
    Tuple[int, int, int, int],
    Tuple[float, float, float, float],
    str
]
# extensions
BIN = ".bin"
GML = ".gml"
LOG = ".log"
MP4 = ".mp4"
PDF = ".pdf"
PKL = ".pkl"
TXT = ".txt"
# paths
load_dotenv()
#
LRGSG_ROOT = Path(os.getenv('LRGSG_ROOT'))
LRGSG_BUILD = Path(os.getenv('LRGSG_BUILD'))
LRGSG_DATA = Path(os.getenv('LRGSG_DATA'))
LRGSG_IPYNB = Path(os.getenv('LRGSG_IPYNB'))
LRGSG_SRC = Path(os.getenv('LRGSG_SRC'))
LRGSG_LIB = Path(os.getenv('LRGSG_LIB'))
LRGSG_LIB_CCORE = Path(os.getenv('LRGSG_LIB_CCORE'))
LRGSG_LIB_CBIN = Path(os.getenv('LRGSG_LIB_CBIN'))
LRGSG_LIB_CONFIG = Path(os.getenv('LRGSG_LIB_CONFIG'))
LRGSG_LIB_GT_PATCHES = Path(os.getenv('LRGSG_LIB_GT_PATCHES'))
LRGSG_LIB_NX_PATCHES = Path(os.getenv('LRGSG_LIB_NX_PATCHES'))
LRGSG_LIB_STOCPROC = Path(os.getenv('LRGSG_LIB_STOCPROC'))
LRGSG_TEST = Path(os.getenv('LRGSG_TEST'))
LRGSG_TOOLS = Path(os.getenv('LRGSG_TOOLS'))
LRGSG_LOG = Path(os.getenv('LRGSG_LOG'))
#
# env_vars = [
#     'LRGSG_ROOT', 'LRGSG_BUILD', 'LRGSG_DATA', 'LRGSG_IPYNB', 'LRGSG_SRC', 'LRGSG_LIB',
#     'LRGSG_LIB_CCORE', 'LRGSG_LIB_CBIN', 'LRGSG_LIB_CONFIG', 'LRGSG_LIB_GT_PATCHES',
#     'LRGSG_LIB_NX_PATCHES', 'LRGSG_LIB_STOCPROC', 'LRGSG_TEST', 'LRGSG_TOOLS', 'LRGSG_LOG'
# ]

# globals().update({var: Path(os.getenv(var)) for var in env_vars})
# #
PATHDATA = Path(LRGSG_DATA)#
#
PATHNLLIB = 'lrgsglib'
#
PATHNDATA = PATHDATA.name
#
PATHNCLDT = 'cluster_data'
PATHNGRPH = 'graph'
PATHNISNG = 'ising'
PATHNLRGS = 'lrgsg'
PATHNPAPR = 'paper'
PATHNPHTR = 'phtra'
PATHNPLOT = 'plot'
PATHNSPEC = 'spect'
PATHNVM = 'voter'
#
PATHPLOT =  PATHDATA / Path(PATHNPLOT)
PATHPAPER = PATHPLOT / Path(PATHNPAPR)
PATHCLDAT = Path(PATHNDATA, PATHNCLDT)
#
# default values
LB_PFLIP = 0.
UB_PFLIP = 1.
#
BSP_VAL = [-1, +1]
BSP_STORE_MODE = ''
BSP_STORE_MODES = [BSP_STORE_MODE, 'persistent', 'sequential']
BSP_STORE_FREQ = 0
BSP_RUN_MODE_C = 'C'
BSP_RUN_MODE_PYTHON = 'Python'
BSP_RUN_MODE = BSP_RUN_MODE_PYTHON
BSP_RUN_MODES_C_LIST = [BSP_RUN_MODE_C, "c"]
BSP_RUN_MODES_C_DICT = {bsp: BSP_RUN_MODE_C for bsp in BSP_RUN_MODES_C_LIST}
BSP_RUN_MODES_PY_LIST = [BSP_RUN_MODE_PYTHON, "py", "Py", "PY", "Python"]
BSP_RUN_MODES_PY_DICT = {bsp: BSP_RUN_MODE_PYTHON for bsp in BSP_RUN_MODES_PY_LIST}
BSP_RUN_MODES = BSP_RUN_MODES_PY_LIST + BSP_RUN_MODES_C_LIST
BSP_RUN_MODES_DICT = {**BSP_RUN_MODES_PY_DICT, **BSP_RUN_MODES_C_DICT}
# Signed Graph default values
SG_PFLIP = 0.0
SG_IMPORT_ON = False
SG_INIT_NW_DICT = False
SG_INIT_WVAL = 1.0
SG_XPRT_MODE = 'pk'
SG_MPRT_MODE = 'pk'
SG_LIST_REPR = ['H']
SG_GRAPH_REPR = 'G'
SG_GRAPHINT_REPR = SG_GRAPH_REPR
SG_GRAPHGEO_REPR = 'H'
SG_ERRMSG_MAXEIGVIDX = """The maximum eigenvalue index is out of bounds."""
SG_ERRMSG_NW_DICT = f"Inheriting class must have attribute 'nwContainer'"
SG_ERRMSG_NFLIP = """The probability of flipping an edge times the 
                             number of edges is < 1, then no edges would be
                             flipped. No flip will be performed."""
SG_ERRMSG_PFLIP = f""" pflip must be between {LB_PFLIP} and {UB_PFLIP}, inclusive."""
SG_ERRMSG_NOCLUST_GTR_NCB = """Requested number of Cluster files is bigger than the one in selected topology."""
SG_ERRMSG_ZEROCLUST = """No clusters were found."""
# erdos renyi default values
WS_PHTABB = "ws"
WS_ONREP = 'G'
WS_STDFN = ""
WS_SGPATH = ""
#
FC_N = 100
# erdos renyi default values
ER_PHTABB = "er"
ER_ONREP = 'G'
ER_STDFN = ""
ER_SGPATH = ""
# 2D Lattice default values
L2D_FBCV = 1.
L2D_PREW = 0.
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
L2D_GEO_SQRSW = 'squared_sw'
L2D_GEO_TRISW = 'triangular_sw'
L2D_GEO_TRI_SHRT = 'tri'
L2D_GEO_SQR_SHRT = 'sqr'
L2D_GEO_HEX_SHRT = 'hex'
L2D_GEO_SQRSW_SHRT = 'sqr_sw'
L2D_GEO_TRISW_SHRT = 'tri_sw'
L2D_P_C_LIST = [0.146, 0.103, 0.065, np.nan, np.nan] 
L2D_WITH_POS = False
# 
L2D_GEO = L2D_GEO_SQR
# 
L2D_GEO_LIST = [L2D_GEO_TRI, L2D_GEO_SQR, L2D_GEO_HEX, L2D_GEO_SQRSW, L2D_GEO_TRISW]
L2D_GEO_SHRT_LIST = [L2D_GEO_TRI_SHRT, L2D_GEO_SQR_SHRT, L2D_GEO_HEX_SHRT, L2D_GEO_SQRSW_SHRT, L2D_GEO_TRISW_SHRT]
L2D_SINGLE_CELL_LIST = ['single', 'singleXERR', 'singleZERR']
L2D_RAND_CELL_LIST = ['rand', 'randXERR', 'randZERR']
#
L2D_P_C_DICT = {g: p for g,p in zip(L2D_GEO_LIST, L2D_P_C_LIST)}
L2D_PATH_DICT = {g: L2D_PHTABB + g for g in L2D_GEO_LIST}
L2D_ZIP_GEO_SHRT = list(zip(L2D_GEO_SHRT_LIST, L2D_GEO_LIST))
L2D_GEO_SHRT_DICT = {s: a for a,s in L2D_ZIP_GEO_SHRT}
L2D_SHRT_GEO_DICT = {a: s for a,s in L2D_ZIP_GEO_SHRT}
#
L2D_WARNMSG_GEO = """The selected geometry of the 2D lattice is not available. Setting it to 'squared' for a 2d regular grid."""
L2D_ERRMSG_GEO = """Invalid side value for hexagonal lattice. In order to implement PBC on hexagonal lattice you need to provide an even value for side1 and side2."""
#
L3D_DIM0 = 16
L3D_DIM = tuple(L3D_DIM0 for _ in range(3))
L3D_PBC = True
L3D_FBCV = 1.
L3D_PDIL = 0.
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
L3D_ONREP = 'G'
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
LRSG_ENTROPY_STEP = 1000
DEFAULT_ENTROPY_LEXPONENT = -3
DEFAULT_ENTROPY_HEXPONENT = 6
DEFAULT_ISING_NOCLUST = 5
DEFAULT_MAX_DIGITS_ROUND_SIGFIG = 18
DEFAULT_NUNMBER_AVERAGES = 100
DEFAULT_SPIKE_THRESHOLD = 0.05
DEFAULT_MAX_THRESHOLD = 2 * DEFAULT_SPIKE_THRESHOLD



COUNT_XERR_PATTERNS = True


op_map = {
    "<":  operator.lt,
    "<=": operator.le,
    ">":  operator.gt,
    ">=": operator.ge,
    "==": operator.eq,
    "!=": operator.ne
}