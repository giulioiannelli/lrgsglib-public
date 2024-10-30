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
out_suffix = args.out_suffix
NoClust = args.NoClust
runlang = args.runlang
remove_files = args.remove_files
workdir = args.workdir
print_chrono = args.print_chrono
if ic.startswith('ground_state'):
    parts = ic.split('_')
    number = int(parts[-1])
    howmany = number+1
else:
    number  = 0
    howmany = 1
#
erDictArgs = dict(n=N, p=p, sgpath=workdir, pflip=pflip, init_nw_dict=True)
isingDictArgs = dict(T=T, ic=ic, runlang=runlang, NoClust=NoClust, rndStr=True, 
                     out_suffix=out_suffix, id_string=in_suffix)
# er = ErdosRenyi(**erDictArgs)
# Fname = f'GC_meanVar_p={pflip:.3g}_{cell}_{in_suffix}.txt'
# pathFname = os.path.join(er.expOut, Fname)
# if not os.path.exists(pathFname):
#     lenList = []
#     for _ in range(navg2):
#         er = ErdosRenyi(**erDictArgs)
#         er.flip_sel_edges(er.nwDict[cell]['G'])
#         # er.flip_random_fract_edges()
#         er.compute_k_eigvV(howmany=number + 1)
#         er.load_eigV_on_graph(which=number, binarize=True)
#         er.make_clustersYN(f'eigV{number}', +1)
#         lenList.append(len(er.gc))
#     meanN, stdN = np.mean(lenList)/er.N, np.std(lenList)
#     np.savetxt(pathFname, np.atleast_2d([meanN, stdN]), fmt='%.3g')
# else:
#     meanN, stdN = np.loadtxt(pathFname)
for _ in range(navg):
    # iter_gen = 0
    # while True:
    er = ErdosRenyi(**erDictArgs)
    er.flip_sel_edges(er.nwDict[cell]['G'])
    # er.flip_random_fract_edges()
    er.compute_k_eigvV(k=howmany)
    er.load_eigV_on_graph(which=number, binarize=True)
    er.make_clustersYN(f'eigV{number}', +1)
        # if  abs(len(er.gc)/er.N - meanN) < stdN:
        #     break
        # elif iter_gen > DEFAULT_MAX_ITER_ER_GC:
        #     raise ValueError(f"Exceeded maximum number of iterations\
        #                       {DEFAULT_MAX_ITER_ER_GC}")
        # iter_gen += 1
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