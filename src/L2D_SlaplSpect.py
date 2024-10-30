from parsers.L2D_SlaplSpect import *
#
args = parser.parse_args()
#
side = args.L
p = args.p
geo = args.geometry
cell = args.cell_type
mode = args.mode
eigmode = args.eigen_mode
navg = args.number_of_averages
period = args.period
bins_count = args.bins_count
howmany = args.howmany
#
#
l = Lattice2D(side1=side, pflip=p, geo=geo, sgpath=args.workDir)
workDir = l.spectpath
#
#
match mode:
    case "eigvec_dist":
        fname_base = '_'.join([f"dist{howmany}", f"{p:.3g}", eigmode])
        fname = workDir+fname_base+f"_{navg}.pkl"
        initial_eig = np.abs(eigV_for_lattice2D_ptch(side=side, pflip=p, geo=geo, 
                                                     mode=eigmode, howmany=howmany))
        bins, bin_centers = create_symmetric_log_bins(np.min(initial_eig), 
                                                     np.max(initial_eig), 
                                                     bins_count)
        bin_counter = [Counter() for i in range(howmany)]
        fname_old = ""
        if not os.path.exists(fname):
            navg_done = 0
            globfname_exists = glob.glob(f"{workDir}{fname_base}*")
            if globfname_exists:
                navg_done = int(os.path.splitext(globfname_exists[0].split('_')[-1])[0])
                fname_old = globfname_exists[0]
            nAvgNeed = navg - navg_done
            for current_period in range((nAvgNeed // period) + bool(nAvgNeed % period)):
                batch_size = min(nAvgNeed - current_period * period, period)
                eig_values = [[] for i in range(howmany)]
                for _ in range(batch_size):
                    eigV = eigV_for_lattice2D_ptch(side=side, pflip=p, geo=geo, 
                                                   mode=eigmode, howmany=howmany)
                    for i in range(howmany):
                        eig_values[i].append(eigV[i])
                eig_values = [np.concatenate(i) for i in eig_values]
                for i in range(howmany):
                    bin_counter[i].update(bin_eigenvalues(eig_values[i], bins, bin_centers))
                fnameNew = f"{workDir}{fname_base}_{navg_done + (current_period + 1) * period}.pkl"
                try:
                    os.rename(fname_old, fnameNew)
                except FileNotFoundError:
                    pass
                with open(fnameNew, "wb") as f:
                    pk.dump(bin_counter, f)
                fname_old = fnameNew

            # At the end, save the final state if needed
            if nAvgNeed % period:
                fnameNew = workDir+fname_base+f"_{navg}.pkl"
                os.rename(fname_old, fnameNew)
                with open(fnameNew, "wb") as f:
                    pk.dump(bin_counter, f)
    case 'eigval_dist':
        fname_old = ''
        fname_base = '_'.join([f"dist_eigval", f"{p:.3g}", cell])
        fname_suffix = f"_{navg}.pkl"
        fname = workDir+fname_base+fname_suffix
        if not os.path.exists(fname):
            navg_done = 0
            globfname_exists = glob.glob(f"{workDir}{fname_base}*")
            if globfname_exists:
                navg_done = int(os.path.splitext(globfname_exists[0].split('_')[-1])[0])
                fname_old = globfname_exists[0]
            nAvgNeed = navg - navg_done
            #
            all_eigv = []
            for current_period in range((nAvgNeed // period) + bool(nAvgNeed % period)):
                batch_size = min(nAvgNeed - current_period * period, period)
                for _ in range(batch_size):
                    l = Lattice2D(side, geo=geo, pflip=p, init_nw_dict=True)
                    l.flip_sel_edges(l.nwDict[cell]['G'])
                    l.compute_full_laplacian_spectrum()
                    all_eigv.extend(l.eigv)
                fnameNew = f"{workDir}{fname_base}_{navg_done + (current_period + 1) * period}.pkl"
                try:
                    os.rename(fname_old, fnameNew)
                except FileNotFoundError:
                    pass
                with open(fnameNew, "wb") as f:
                    pk.dump(all_eigv, f)
                fname_old = fnameNew

            # At the end, save the final state if needed
            if nAvgNeed % period:
                fnameNew = workDir+fname_base+f"_{navg}.pkl"
                os.rename(fname_old, fnameNew)
                with open(fnameNew, "wb") as f:
                    pk.dump(all_eigv, f)
