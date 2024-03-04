from Lattice2D_phasetr_Parser import *
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

Pinf=np.zeros(args.number_of_averages)
Pinf2=np.zeros(args.number_of_averages)
Fluct=np.zeros(args.number_of_averages)
Fluct2=np.zeros(args.number_of_averages)

lattice = Lattice2D(args.L, pflip=args.p, geometry=args.geometry)
file1 = f'{lattice.phtrapath}p_{lattice.pflip}_{args.cell_type}_{args.number_of_averages-1}'
print(file1)

if os.path.exists(file1):
    os.exit(0)

for cont, avg in enumerate(range(args.number_of_averages)):
    lattice = Lattice2D(args.L, pflip=args.p, geometry=args.geometry)
    lattice.flip_sel_edges(geometry_func(lattice))

    file1 = f'{lattice.phtrapath}p_{lattice.pflip}_{args.cell_type}_{avg}'

    lattice.compute_k_eigvV()
    eigV = lattice.eigV[0]
    lattice.calc_fluct_Pinf()

    Fluct[cont]=lattice.eigV_fluct
    Fluct2[cont]=lattice.eigV_fluct**2
    
    Pinf[cont]=lattice.Pinf
    Pinf2[cont]=lattice.Pinf**2
        
    x=[np.sum(Pinf)/(1+avg),np.sum(Fluct)/(1+avg), np.sum(Fluct2)/(1+avg), np.var(Fluct[Fluct!=0]), np.sum(Pinf2)/(1+avg),lattice.pflip,int(avg+1)]
    np.savetxt(file1, np.atleast_2d(x), fmt='%.7g', delimiter=',')
