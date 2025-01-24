from lrgsglib.core import *

def initialize_l2d_dict_args(args):
    return dict(side1=args.L, geo=args.geometry, sgpathn=args.workdir, 
                pflip=args.p, init_nw_dict=True)

def initialize_ising_dict_args(args, out_suffix):
    return dict(T=args.T, ic=args.init_cond, runlang=args.runlang, 
                NoClust=args.NoClust, rndStr=True,
                out_suffix=out_suffix)

def get_out_suffix(args, ic_gs, number):
    return args.out_suffix or "gs"+str(number)+args.cell_type if ic_gs else args.init_cond

def run_simulation(args, l2dDictArgs, isingDictArgs, number, remove_files, min_gc_cond=True):
    MAX_COUNT = 1e5
    for _ in range(args.number_of_averages):
        if args.runlang == 'C4':
            lengc = 0
            count = 0
            while lengc < 0.33 or count < MAX_COUNT:
                seedsave = np.random.randint(0, 2**32-1)
                lattice = Lattice2D(seed=seedsave, **l2dDictArgs)
                lattice.flip_sel_edges(lattice.nwDict['randXERR']['G'])
                lattice.compute_k_eigvV(number+1)
                lattice.load_eigV_on_graph(number, binarize=True)
                lattice.make_clustersYN(f"eigV{number}", val=-1)
                lengc = len(lattice.gc)/lattice.N
                count += 1
            
        else:
            lattice = Lattice2D(**l2dDictArgs)
            lattice.flip_sel_edges(lattice.nwDict[args.cell_type]['G'])
            lattice.compute_k_eigvV(k=number+1)
        #
        isdy = IsingDynamics(lattice, **isingDictArgs)
        isdy.init_ising_dynamics()
        lattice.export_eigV(number, exName=isdy.id_string_isingdyn)
        lattice.export_edgel_bin(exName=isdy.id_string_isingdyn)
        isdy.export_ising_clust(which=number, val=-1)
        isdy.run(verbose=False, thrmSTEP=args.thrmsteps)
        if remove_files:
            isdy.remove_run_c_files(remove_stderr=True)
            lattice.remove_edgl_file()
            if args.runlang == 'C4':
                lattice.remove_eigV_file()
