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
DEFAULT_INITCON = "ground_state_0"
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
    "Nclust",
    default=5,
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
    "-nA",
    "-navg",
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
parser.add_argument(
    "--initial_cond",
    default=DEFAULT_INITCON,
    type=str,
)

args = parser.parse_args()
#
import_on = False
if os.path.exists(
    DEFAULT_GRAPH_OUTDIR
    + f"N={args.L**2:d}/"
    + args.graph_filename
    + f"_p={args.p:.3g}.pickle"
):
    import_on = True
#
sqlatt = Lattice2D(
    side1=args.L,
    geometry="squared",
    import_on=import_on,
    pflip=args.p,
    stdFname=args.graph_filename + f"_p={args.p:.3g}",
    expathc= f"N={args.L**2:}/"
)
if not import_on:
    sqlatt.flip_random_fract_edges()
    sqlatt.export_graph_pickle()
    sqlatt.export_adj_bin()

ising_dyn = IsingDynamics(
    system=sqlatt, 
    T=args.T, 
    IsingIC=args.initial_cond,
    MODE_RUN="C1", 
    NoClust = 5
)

ising_dyn.init_ising_dynamics()
if not os.path.exists(
    DEFAULT_GRAPH_OUTDIR
    + "s_"
    + args.graph_filename
    + f"_p={args.p:.3g}.pickle"
):
    ising_dyn.export_s_init()
if args.p < 0.103:
    if not os.path.exists(
        DEFAULT_GRAPH_OUTDIR
        + "cl1_"
        + args.graph_filename
        + f"_p={args.p:.3g}.pickle"
    ):
        ising_dyn.find_ising_clusters()
        ising_dyn.mapping_nodes_to_clusters()
        ising_dyn.export_ising_clust(howmany=1)
for _ in range(args.number_of_averages):
    ising_dyn.run(out_suffix=args.out_suffix, tqdm_on=False)

