from parsers.Lattice2D_SlapSpect_Parser import *
#
geo = 'squared'
nAvg = 500
bins_count = 500
period = 20
workDir = "data/plot/dev_tool_Lattice2D/"
mode = 'scipy'
side = 18
pflip = 0.01
howmany = 1
fnameBase = f"dist_{side}_{pflip:.3g}_{mode}"
fname = f"{workDir}{fnameBase}_{nAvg}.pkl"
#
os.makedirs(workDir, exist_ok=True)
#
initial_eig = np.abs(eigV_for_lattice2D_ptch(side, pflip=pflip, geo=geo, 
                                             mode=mode, howmany=howmany))
bins, bin_centers = create_symmetric_log_bins(np.min(initial_eig), 
                                              np.max(initial_eig), 
                                              bins_count)
bin_counter = [Counter() for i in range(howmany)]
fnameOld = ""
if not os.path.exists(fname):
    nAvgDone = 0
    fnameExists = glob.glob(f"{workDir}{fnameBase}*")
    if fnameExists:
        nAvgDone = int(os.path.splitext(fnameExists[0].split('_')[-1])[0])
        fnameOld = fnameExists[0]
    nAvgNeed = nAvg - nAvgDone
    for current_period in range((nAvgNeed // period) + bool(nAvgNeed % period)):
        batch_size = min(nAvgNeed - current_period * period, period)
        eig_values = [[] for i in range(howmany)]
        for _ in range(batch_size):
            eigV = eigV_for_lattice2D_ptch(side, pflip=pflip, geo=geo, 
                                                mode=mode, howmany=howmany)
            for i in range(howmany):
                eig_values[i].append(eigV[i])
        eig_values = [np.concatenate(i) for i in eig_values]
        for i in range(howmany):
            bin_counter[i].update(bin_eigenvalues(eig_values[i], bins, bin_centers))
        fnameNew = f"{workDir}{fnameBase}_{nAvgDone + (current_period + 1) * period}.pkl"
        try:
            os.rename(fnameOld, fnameNew)
        except FileNotFoundError:
            pass
        with open(fnameNew, "wb") as f:
            pk.dump(dict(bin_counter), f)
        fnameOld = fnameNew

    # At the end, save the final state if needed
    if nAvgNeed % period:
        fnameNew = f"{workDir}{fnameBase}_{nAvg}.pkl"
        os.rename(fnameOld, fnameNew)
        with open(fnameNew, "wb") as f:
            pk.dump(dict(bin_counter), f)
