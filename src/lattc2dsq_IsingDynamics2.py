from LRGSG_package.LRGSG import *

#
description = """
    Compute average magnetization of a square lattice with prescribed adjacency 
    matrix.
"""
DEFAULT_NUNMBER_AVERAGES = 100
DEFAULT_GRAPH_OUTDIR = "src/LRGSG_package/dump/"
DEFAULT_GRAPH_NAME = "sqLattice"
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
    type=int,
)

args = parser.parse_args()
#
import_on = False
if os.path.exists(DEFAULT_GRAPH_OUTDIR + args.gn + f"_p={args.p:.3g}.pickle"):
    import_on = True
sqlatt = Lattice2D(
    side1 = args.L,
    geometry = 'squared',
    import_on = import_on,
    pflip = args.p,
    stdFname = args.gn + f"_p={args.p:.3g}"
)
if not import_on:
    sqlatt.flip_random_fract_edges()
    sqlatt.export_graph_pickle()
    sqlatt.export_adj_bin()

ising_dyn = IsingDynamics(
    system = sqlatt,
    T = args.T, 
    IsingIC = 'ground_state'
)
ising_dyn.init_ising_dynamics()
ising_dyn.run(MODE="C")


ising_dyn.id_string_isingdyn
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