from parsers.Lattice2D_TransCluster_Parser import *
#
args = parser.parse_args()
if args.cell_type in ['rand', 'randZERR', 'randXERR']:
    geometry_func = lambda lattice: lattice.neg_weights_dict['randZERR']['G']
elif args.cell_type == 'cross':
    geometry_func = lambda lattice: lattice.neg_weights_dict['randXERR']['G']
elif args.cell_type == 'single':
    geometry_func = lambda lattice: lattice.neg_weights_dict['rand']['G']
else:
    raise ValueError("Invalid cell specified")
#
lattice = Lattice2D(args.L, pflip=args.p, geo=args.geometry, 
                    init_weight_dict=False, with_positions=False)
mpath = {'pCluster': lattice.lrgsgpath,
       'ordParam': lattice.phtrapath}
#
if args.mode == 'pCluster':
    merged_dict = Counter()
    extout = ePKL
elif args.mode == 'ordParam':
    Pinf = np.zeros(args.number_of_averages)
    Pinf2 = np.zeros(args.number_of_averages)
    Fluct = np.zeros(args.number_of_averages)
    Fluct2 = np.zeros(args.number_of_averages)
    extout = eTXT
#
def file_path_maker(mpath, ppath = args.p, 
                    napath = args.number_of_averages, 
                    spath = args.out_suffix,
                    ctpath = args.cell_type,
                    extout = extout):
    return f'{mpath}{args.mode}_p={ppath:.3g}_{ctpath}_na={napath}{spath}{extout}'
#

#
if args.mode == 'pCluster':
    filename = file_path_maker(mpath[args.mode])
    if os.path.exists(filename):
        exit(f"File {os.path.split(filename)[1]} already exists.")
    for avg in range(args.number_of_averages):
        lattice = Lattice2D(args.L, pflip=args.p, geo=args.geometry)
        lattice.flip_sel_edges(geometry_func(lattice))
        #
        dist_dict = lattice.cluster_distribution_list()
        merged_dict += Counter(dist_dict)
        #
        if avg % 100 == 0:
            try:
                filenameold = file_path_maker(mpath[args.mode], napath=avg)
                os.remove(filenameold)
            except OSError:
                pass
            filename = file_path_maker(mpath[args.mode], napath=avg+100)
            with open(filename, 'wb') as file:
                pickle.dump(merged_dict, file)
elif args.mode == 'ordParam':
    for cont, avg in enumerate(range(args.number_of_averages)):
        lattice = Lattice2D(args.L, pflip=args.p, geo=args.geometry)
        lattice.flip_sel_edges(geometry_func(lattice))
        #
        lattice.compute_k_eigvV()
        lattice.calc_fluct_Pinf()
        #
        Fluct[cont]=lattice.eigV_fluct
        Fluct2[cont]=lattice.eigV_fluct**2
        Pinf[cont]=lattice.Pinf
        Pinf2[cont]=lattice.Pinf**2
        #
        data=[lattice.pflip,
              lattice.Ne_n,
              avg+1,
              np.sum(Pinf), 
              np.sum(Pinf2), 
              np.sum(Fluct), 
              np.sum(Fluct2), 
              np.var(Fluct[Fluct!=0])]
        #
        if avg % 100 == 0:
            try:
                filenameold = file_path_maker(mpath[args.mode], napath=avg)
                os.remove(filenameold)
            except OSError:
                pass
            filename = file_path_maker(mpath[args.mode], napath=avg+100)
            with open(filename, 'wb') as file:
                np.savetxt(file, np.atleast_2d(data), fmt='%.7g')