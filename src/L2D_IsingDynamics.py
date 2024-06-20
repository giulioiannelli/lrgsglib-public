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
in_suffix = args.in_suffix
out_suffix = args.out_suffix
NoClust = args.NoClust
runlang = args.runlang
remove_files = args.remove_files
if ic.startswith('ground_state'):
    parts = ic.split('_')
    number = int(parts[-1])
    howmany = number+1
else:
    howmany = 1
#
for na in range(navg):
    l = Lattice2D(side1=side, geo=geo, pflip=p, init_nw_dict=True)
    l.flip_sel_edges(l.nwDict[cell]['G'])
    l.compute_k_eigvV(howmany=howmany)
    isdy = IsingDynamics(l, T, ic=ic, runlang=runlang, NoClust=NoClust, 
                            rndStr=True, out_suffix=out_suffix, 
                            id_string=in_suffix)
    isdy.init_ising_dynamics()
    l.export_edgel_bin(exName=isdy.id_string_isingdyn)
    isdy.export_ising_clust()
    isdy.run(verbose=False)
    if remove_files:
        isdy.remove_run_c_files(remove_stderr=True)
        l.remove_edgl_file()
