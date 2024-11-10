from parsers.ER_IsingDynamics import *
#
args = parser.parse_args()
#
N = args.N
p = args.p
pflip = args.pflip
T = args.T
# geo = args.geometry
thrmsteps = args.thrmsteps
cell = args.cell_type
ic = args.init_cond
navg = args.number_of_averages
navg2 = args.number_of_averages2
in_suffix = args.in_suffix
NoClust = args.NoClust
runlang = args.runlang
remove_files = args.remove_files
workdir = args.workdir
print_chrono = args.print_chrono
ic_gs = ic.startswith('ground_state')
number = int(ic.split('_')[-1]) if ic_gs else 0
out_suffix = args.out_suffix or "gs"+str(number) if ic_gs else ic

#
erDictArgs = dict(n=N, p=p, sgpath=workdir, pflip=pflip, init_nw_dict=True)
isingDictArgs = dict(T=T, ic=ic, runlang=runlang, NoClust=NoClust, rndStr=True, 
                     out_suffix=out_suffix, id_string=in_suffix)
for _ in range(navg):
    er = ErdosRenyi(**erDictArgs)
    er.flip_sel_edges(er.nwDict[cell]['G'])
    er.compute_k_eigvV(k=number+1)
    er.load_eigV_on_graph(which=number, binarize=True)
    er.make_clustersYN(f'eigV{number}', +1)
    #
    isdy = IsingDynamics(er, **isingDictArgs)
    isdy.init_ising_dynamics()
    er.export_edgel_bin(exName=isdy.id_string_isingdyn)
    isdy.export_ising_clust()
    isdy.run(verbose=False, thrmSTEP=thrmsteps)
    if remove_files:
        isdy.remove_run_c_files(remove_stderr=True)
        er.remove_edgl_file()
if print_chrono:
    Chronometer.print_all_chronometers()

print("Done!")