from typing import Any, Optional, Union, List
from scipy.sparse import csr_array
from networkx.classes.graph import Graph
# extensions
ePDF = ".pdf"
eTXT = ".txt"
eBIN = ".bin"
eMP4 = ".mp4"
ePKL = ".pkl"
# paths
PATH_ROOTF = "LRG-Signed"
# default values
DEFLBpflip = .0
DEFUBpflip = 1.


DIR_SRC_DEFAULT = "src/"
DIR_PCK_DEFAULT = "LRGSG_package/"
DIR_PLT_DEFAULT = "plot/"
DIR_DAT_DEFAULT = "data/"

DIR_GRAPH_DEFAULT = "graphs/"
DIR_ISING_DEFAULT = "ising/"
DIR_VOTER_DEFAULT = "voter/"
DIR_LRGSG_DEFAULT = "lrgsg/"
DIR_PHTRA_DEFAULT = "phtra/"


DEFErdosReny_pthabb = "er/"

DEFLattice2D_side1 = 32
DEFLattice2D_side2 = 0
DEFLattice2D_pbc = True
DEFLattice2D_fbcv = 1.
DEFLattice2D_stdFn = ""
DEFLattice2D_sgpath = ""
DEFLattice2D_pthabb = 'l2d_'
DEFLattice2D_geoTri = 'triangular'
DEFLattice2D_geoTri_abb = 'tri'
DEFLattice2D_geoSqr = 'squared'
DEFLattice2D_geoSqr_abb = 'sqr'
DEFLattice2D_geoHex = 'hexagonal'
DEFLattice2D_geoHex_abb = 'hex'
DEFLattice2D_p_clist = [0.146, 0.103, 0.065] 
#
DEFLattice2D_geo = DEFLattice2D_geoSqr
#
DEFLattice2D_geoabblist = [DEFLattice2D_geoTri_abb, 
                        DEFLattice2D_geoSqr_abb, 
                        DEFLattice2D_geoHex_abb]
DEFLattice2D_geolist = [DEFLattice2D_geoTri, 
                        DEFLattice2D_geoSqr, 
                        DEFLattice2D_geoHex]
#

DEFLattice2D_p_cdict = {g: p for g,p in 
                        zip(DEFLattice2D_geolist, DEFLattice2D_p_clist)}
DEFLattice2D_pthdict = {g: f"{DEFLattice2D_pthabb}{g}/" for a,g in 
                        zip(DEFLattice2D_geoabblist, DEFLattice2D_geolist)}
DEFLattice2D_geodict = {s: a for a,s in 
                        zip(DEFLattice2D_geoabblist, DEFLattice2D_geolist)}
DEFLattice2D_geodictabb = {a: s for a,s in 
                        zip(DEFLattice2D_geoabblist, DEFLattice2D_geolist)}

DEFLattice2D_geowarnmsg = """The selected geometry of the 2D lattice is not available. Setting it to 'squared' for a 2d regular grid."""
DEFLattice2D_geoerrmsg = """Invalid side value for hexagonal lattice. In order to implement PBC on hexagonal lattice you need to provide an even value for side1 and side2."""
DEFSignedGraph_pflipverr = f""" pflip must be between {DEFLBpflip} and {DEFUBpflip}, inclusive."""









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
