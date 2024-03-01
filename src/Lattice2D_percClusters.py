from Lattice2D_percClusters_Parser import *
import faulthandler

faulthandler.enable()
#
args = parser.parse_args()
#
merged_dict = Counter()
#
lattice = Lattice2D(args.L, pflip=args.p, geometry=args.geometry)
def file_path_maker(mpath = lattice.lrgsgpath, ppath= args.p, 
                    napath = args.number_of_averages, 
                    spath = args.out_suffix):
    return f'{mpath}p={ppath:.3g}_na={napath}_{spath}.pkl'
#
filename = file_path_maker()
if os.path.exists(filename):
    os.exit(0)
#
for avg in range(args.number_of_averages):
    lattice = Lattice2D(args.L, pflip=args.p, geometry=args.geometry)
    lattice.flip_random_fract_edges()
    dist_dict = lattice.cluster_distribution_list()
    merged_dict += Counter(dist_dict)
    try:
        filenameold = file_path_maker(napath=avg)
        os.remove(filenameold)
    except OSError:
        pass
    filename = file_path_maker(napath=avg+1)
    with open(filename, 'wb') as file:
        pickle.dump(merged_dict, file)