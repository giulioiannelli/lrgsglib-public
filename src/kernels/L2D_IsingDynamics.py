from lrgsglib.core import *

def initialize_l2d_dict_args(args):
    return dict(side1=args.L, geo=args.geometry, sgpathn=args.workdir, 
                pflip=args.p, init_nw_dict=True)

def initialize_ising_dict_args(args, out_suffix, NoClust):
    return dict(T=args.T, ic=args.init_cond, runlang=args.runlang, 
                NoClust=NoClust, rndStr=True, freq=args.freq,
                out_suffix=out_suffix, thrmSTEP=args.thrmsteps)

def get_out_suffix(args):
    return join_non_empty('_', args.init_cond, args.cell_type, args.out_suffix)

def run_simulation(args):
    ic_gs = args.init_cond.startswith('ground_state') or args.init_cond.startswith('gs')
    number = int(args.init_cond.split('_')[-1]) if ic_gs else 0
    val = ConditionalPartitioning(args.val)
    for _ in range(args.number_of_averages):
        lattice = Lattice2D(**initialize_l2d_dict_args(args))
        lattice.flip_sel_edges(lattice.nwDict[args.cell_type]['G'])
        lattice.compute_k_eigvV(number+1)
        lattice.load_eigV_on_graph(number, binarize=True)
        lattice.make_clustersYN(f'eigV{number}', val=val)
        #
        NoClust = lattice.handle_no_clust(NoClust=args.NoClust)
        isingDictArgs = initialize_ising_dict_args(args, get_out_suffix(args), NoClust)
        #
        if NoClust is None:
            continue
        #
        isdy = IsingDynamics(lattice, **isingDictArgs)
        isdy.init_ising_dynamics(exName=isdy.id_string_isingdyn)
        match args.runlang:
            case 'C4'|'C1':
                lattice.export_ising_clust(exName=isdy.id_string_isingdyn, 
                                        NoClust=NoClust)
                if args.runlang == 'C4':
                    lattice.export_eigV(number, exName=isdy.id_string_isingdyn, verbose=args.verbose)
        isdy.run(verbose=args.verbose)
        if args.remove_files:
            isdy.remove_run_c_files(remove_stderr=True)
            lattice.remove_exported_files()