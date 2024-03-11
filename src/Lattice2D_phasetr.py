from Lattice2D_phasetr_Parser import *
#
args = parser.parse_args()
if args.cell_type in ['square', 'triangle', 'hexagon']:
    geometry_func = lambda lattice: lattice.neg_weights_dict.nwdZERR
elif args.cell_type == 'single':
    geometry_func = lambda lattice: lattice.neg_weights_dict.nwd_reflip
elif args.cell_type == 'cross':
    geometry_func = lambda lattice: lattice.neg_weights_dict.nwdXERR
else:
    raise ValueError("Invalid cell specified")
#
Pinf = np.zeros(args.number_of_averages)
Pinf2 = np.zeros(args.number_of_averages)
Fluct = np.zeros(args.number_of_averages)
Fluct2 = np.zeros(args.number_of_averages)

lattice = Lattice2D(args.L, pflip=args.p, geo=args.geometry)
def file_path_maker(mpath = lattice.phtrapath, ppath= args.p, 
                    napath = args.number_of_averages, 
                    spath = args.out_suffix,
                    ctpath = args.cell_type):
    return f'{mpath}p={ppath:.3g}_{ctpath}_na={napath}_{spath}.txt'


if os.path.exists(file_path_maker()):
    os.exit(0)


for cont, avg in enumerate(range(args.number_of_averages)):
    lattice = Lattice2D(args.L, pflip=args.p, geo=args.geometry)
    try:
        filenameold = file_path_maker(napath=avg)
        os.remove(filenameold)
    except OSError:
        pass
    lattice.flip_sel_edges(geometry_func(lattice))


    lattice.compute_k_eigvV()
    eigV = lattice.eigV[0]
    lattice.calc_fluct_Pinf()

    Fluct[cont]=lattice.eigV_fluct
    Fluct2[cont]=lattice.eigV_fluct**2
    
    Pinf[cont]=lattice.Pinf
    Pinf2[cont]=lattice.Pinf**2
        
    x=[lattice.pflip, np.sum(Pinf)/(1+avg), np.sum(Pinf2)/(1+avg), np.sum(Fluct)/(1+avg), np.sum(Fluct2)/(1+avg), 
       np.var(Fluct[Fluct!=0]),int(avg+1)]
    
    filename = file_path_maker(napath=avg+1)
    with open(filename, 'wb') as file:
        np.savetxt(file, np.atleast_2d(x), fmt='%.7g')
