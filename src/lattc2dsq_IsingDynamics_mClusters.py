from LRGSG_package.LRGSG import *
from LRGSG_package.LRGSG_utils import move_to_rootf
from lattc2dsq_IsingDynamics_mClusters_Parser import *
#
# move_to_rootf()
# per il programma C salva <m>, <m^2>
# ogni volta sovrascrivi
#
parser = argparse.ArgumentParser(description=description)
parser.add_argument(
    "L",
    help=HELP_L,
    type=int,
)
parser.add_argument(
    "T",
    help=HELP_T,
    type=float,
)
parser.add_argument(
    "p",
    help=HELP_p,
    type=float,
)
parser.add_argument(
    "-nC",
    "--Nclust",
    default=DEFAULT_ISING_NOCLUST,
    help=HELP_Nclust,
    type=int,
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
    "-iC",
    "--initial_cond",
    default=DEFAULT_INITCON,
    type=str,
)
args = parser.parse_args()
#
import_on = False
grphNdir = f"{DEFAULT_GRAPH_OUTDIR}N={args.L**2:d}/"
grphpth = f"{grphNdir}{args.graph_filename}_p={args.p:.3g}.pickle"
if os.path.exists(grphpth):
    import_on = True
#
stdFnameSFFX = "" if args.graph_filename == DEFAULT_GRAPH_NAME else args.graph_filename
sqlatt = Lattice2D(
    side1=args.L, 
    geometry="squared", 
    pflip=args.p, 
    import_on=import_on,
    stdFnameSFFX=stdFnameSFFX
)
if not import_on:
    sqlatt.flip_random_fract_edges()
    sqlatt.export_graph()
    sqlatt.export_adj_bin()
#
ising_dyn = IsingDynamics(
    system=sqlatt,
    T=args.T,
    IsingIC=args.initial_cond,
    MODE_RUN="C1",
    NoClust=args.Nclust,
)
ising_dyn.init_ising_dynamics()
ising_dyn.export_s_init()

pathcl = lambda i: f"{ising_dyn.system.isingpath}cl{i}_{ising_dyn.system.stdFname}.bin"
pathcl_list = [pathcl(i) for i in range(ising_dyn.NoClust)]
cluster_idxs = [os.path.exists(pathcl(i)) for i in range(ising_dyn.NoClust)]

if not all(cluster_idxs):
    ising_dyn.find_ising_clusters()
    ising_dyn.export_ising_clust()

for _ in range(args.number_of_averages):
    ising_dyn.run(tqdm_on=False)
