from Lattice2D_pClusters_Parser import *
#
args = parser.parse_args()
if args.cell_type == 'square':
    geometry_func = lambda lattice: lattice.neg_weights_dict.NEG_WEIGHTS_DICT_H_PSQUARE
elif args.cell_type == 'triangle':
    geometry_func = lambda lattice: lattice.neg_weights_dict.NEG_WEIGHTS_DICT_H_PTRIA
elif args.cell_type == 'hexagon':
    geometry_func = lambda lattice: lattice.neg_weights_dict.NEG_WEIGHTS_DICT_H_PHEXA
elif args.cell_type == 'single':
    geometry_func = lambda lattice: lattice.neg_weights_dict.NEG_WEIGHTS_DICT_H_PFLIP
elif args.cell_type == 'cross':
    geometry_func = lambda lattice: lattice.neg_weights_dict.NEG_WEIGHTS_DICT_H_PCROSS
else:
    raise ValueError("Invalid cell specified")
#
merged_dict = Counter()
#
lattice = Lattice2D(args.L, pflip=args.p, geo=args.geometry)
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
    lattice = Lattice2D(args.L, pflip=args.p, geo=args.geometry)
    lattice.flip_sel_edges(geometry_func(lattice))
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