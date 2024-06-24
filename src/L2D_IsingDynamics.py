from parsers.L2D_IsingDynamics import *
#
args = parser.parse_args()
#
side = args.L
p = args.p
T = args.T
geo = args.geometry
cell = args.cell_type
ic = args.init_cond
navg = args.number_of_averages
navg2 = 500
in_suffix = args.in_suffix
out_suffix = args.out_suffix
NoClust = args.NoClust
runlang = args.runlang
remove_files = args.remove_files
workdir = args.workdir
if ic.startswith('ground_state'):
    parts = ic.split('_')
    number = int(parts[-1])
    howmany = number+1
else:
    number  = 0
    howmany = 1
#
l2dDictArgs = dict(side1=side, geo=geo, sgpath=workdir, pflip=p, 
                   init_nw_dict=True)
isingDictArgs = dict(T=T, ic=ic, runlang=runlang, NoClust=NoClust, rndStr=True, 
                     out_suffix=out_suffix, id_string=in_suffix)
l = Lattice2D(**l2dDictArgs)
fName = os.path.join(l.expOutdir, f'GC_meanVar_p={p:.3g}_{cell}_{in_suffix}.txt')
if not os.path.exists(fName):
    lenList = []
    for _ in range(navg2):
        l = Lattice2D(**l2dDictArgs)
        l.flip_sel_edges(l.nwDict[cell]['G'])
        l.compute_k_eigvV(howmany=number + 1)
        l.load_eigV_on_g(which=number, binarize=True)
        l.make_clustersYN(f'eigV{number}', +1)
        lenList.append(len(l.gc))
    meanN, stdN = np.mean(lenList)/l.N, np.std(lenList)/l.N
    np.savetxt(fName, [meanN, stdN], fmt='%.3g')
else:
    meanN, stdN = np.loadtxt(fName)

for _ in range(navg):
    while True:
        l = Lattice2D(**l2dDictArgs)
        l.flip_sel_edges(l.nwDict[cell]['G'])
        l.compute_k_eigvV(howmany=howmany)
        l.load_eigV_on_g(which=number, binarize=True)
        l.make_clustersYN(f'eigV{number}', +1)
        if  abs(len(l.gc)/l.N - meanN) < stdN:
            break
    isdy = IsingDynamics(l, **isingDictArgs)
    isdy.init_ising_dynamics()
    l.export_edgel_bin(exName=isdy.id_string_isingdyn)
    isdy.export_ising_clust()
    isdy.run(verbose=False)
    if remove_files:
        isdy.remove_run_c_files(remove_stderr=True)
        l.remove_edgl_file()
