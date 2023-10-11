from LRGSG_package.LRGSG import *

#
description = """
    Compute average magnetization of a square lattice with prescribed adjacency 
    matrix.
"""
DEFAULT_NUNMBER_AVERAGES = 100
DEFAULT_DATA_OUTDIR = "data/l2d_sq_ising/"
DEFAULT_GRAPH_OUTDIR = DEFAULT_DATA_OUTDIR + "graphs/"
DEFAULT_GRAPH_NAME = "sqLattice"
DEFAULT_OUT_SUFFIX = 0
HELP_L = f"""
    (int) L of the square lattice
"""
HELP_p = f"""
    (float) fraction of edges to flip
"""
HELP_T = f"""
    (float) temperature of the system
"""

DEFAULT_nA = f"| default={DEFAULT_NUNMBER_AVERAGES}"
HELP_nA = f"""
    (float) temperature of the system {DEFAULT_nA:->10}
"""

DEFAULT_gn = f"| default={DEFAULT_GRAPH_NAME}"
HELP_gn = f"""
    (str) file name of the graph container {DEFAULT_gn:->10}
"""

DEFAULT_o = f"| default={DEFAULT_OUT_SUFFIX}"
HELP_o = f"""
    (str) file name of the graph container {DEFAULT_o:->10}
"""

parser = argparse.ArgumentParser(description=description)
parser.add_argument(
    "L",
    help=HELP_L,
    type=int,
)
parser.add_argument(
    "p",
    help=HELP_p,
    type=float,
)
parser.add_argument(
    "T",
    help=HELP_T,
    type=float,
)
parser.add_argument(
    "-nA", "-navg",
    "--number_of_averages",
    default=DEFAULT_NUNMBER_AVERAGES,
    help=HELP_nA,
    type=int,
)
parser.add_argument(
    "-gn",
    "--graph_filename",
    default=DEFAULT_GRAPH_NAME,
    help=HELP_gn,
    type=str,
)
parser.add_argument(
    "-o",
    "--out_suffix",
    default=DEFAULT_OUT_SUFFIX,
    help=HELP_o,
    type=str,
)

args = parser.parse_args()
#
import_on = False
if os.path.exists(DEFAULT_GRAPH_OUTDIR + args.graph_filename + f"_p={args.p:.3g}.pickle"):
    import_on = True
sqlatt = Lattice2D(
    side1 = args.L,
    geometry = 'squared',
    import_on = import_on,
    pflip = args.p,
    stdFname = args.graph_filename + f"_p={args.p:.3g}"
)
if not import_on:
    sqlatt.flip_random_fract_edges()
    sqlatt.export_graph_pickle()
    sqlatt.export_adj_bin()

ising_dyn = IsingDynamics(
    system = sqlatt,
    T = args.T, 
    IsingIC = 'ground_state',
    MODE_RUN = "C"
)

ising_dyn.init_ising_dynamics()
if not os.path.exists(DEFAULT_GRAPH_OUTDIR + "s_" + args.graph_filename + f"_p={args.p:.3g}.pickle"):
    ising_dyn.export_s_init()
if not os.path.exists(DEFAULT_GRAPH_OUTDIR + "cl1_" + args.graph_filename + f"_p={args.p:.3g}.pickle"):
    ising_dyn.find_ising_clusters()
    ising_dyn.mapping_nodes_to_clusters()
    ising_dyn.export_ising_clust(howmany=1)
for _ in range(args.number_of_averages):
    ising_dyn.run(out_suffix = args.out_suffix, tqdm_on = False)

"""
lrgsg = SignedLaplacianAnalysis(#
    system = sqlatt,
    initCond = 'all_1',
    t_steps = 10,
    no_obs = 200
)

sqlatt.flip_random_fract_edges()
ising_dyn = IsingDynamics(
    system = sqlatt,
    T = 0.2, 
    IsingIC = 'ground_state'
)
if not sqlatt.import_on:
    sqlatt.export_graph_pickle()
    sqlatt.export_adj_bin()
"""