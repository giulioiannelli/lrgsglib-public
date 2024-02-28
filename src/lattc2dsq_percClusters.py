from lattc2dsq_percClusters_Parser import *
#
args = parser.parse_args()
#
merged_dict = Counter()
for avg in range(args.number_of_averages):
    lattice = Lattice2D(args.L, pflip=args.p)
    lattice.flip_random_fract_edges()
    dist_dict = lattice.cluster_distribution_list()
    merged_dict += Counter(dist_dict)
    try:
        filenameold = f'{lattice.lrgsgpath}p={args.p:.3g}_na={avg-1}_{args.out_suffix}.pickle'
        os.remove(filenameold)
    except OSError:
        pass
    filename = f'{lattice.lrgsgpath}p={args.p:.3g}_na={avg}_{args.out_suffix}.pickle'
    with open(filename, 'wb') as file:
        pickle.dump(list(merged_dict.values()), file)