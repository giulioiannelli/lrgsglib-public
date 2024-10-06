from parsers.L3D_TransCluster import *
#
args = parser.parse_args()
#
side = args.L
p = args.p
pdil = args.pdil
geo = args.geometry
cell = args.cell_type
mode = args.mode
navg = args.number_of_averages
sfreq = args.save_frequency if args.save_frequency else navg // 20
outsx = args.out_suffix
typf = args.float_type
#
if mode == 'pCluster':
    merged_dict = Counter()
    extout = PKL
elif mode == 'ordParam':
    extout = TXT
else:
    raise ValueError("Invalid mode specified")
#
def file_path_maker(mpath, mode=mode, ppath = p, napath = navg, spath = outsx,
                    ctpath = cell, extout = extout, pdil = pdil):
    pdilStr = f"pdil={pdil:.3g}" if pdil != 0. else ""
    strName = '_'.join(filter(None, [mode, f"p={ppath:.3g}", ctpath, 
                     f"na={napath}", pdilStr, spath]))
    return os.path.join(mpath, strName) + extout
#
lattice = Lattice3D(dim=tuple(side for _ in range(3)), pflip=p, geo=geo)
#
def geometry_func_rand(lattice: Lattice3D):
    return lattice.nwDict[cell]['G']

def geometry_func_ball(lattice: Lattice3D, radius):
    return lattice.nwDict.get_links_rball(radius)
#
match cell:
    case 'rand'|'randZERR'|'randXERR':
        geometry_func = geometry_func_rand
    case s if s.startswith('ball'):
        radius = get_first_int_in_str(cell)
        geometry_func = lambda lattice: geometry_func_ball(lattice, radius)
    case _:
        raise ValueError("Invalid cell specified")
mpath = {'pCluster': lattice.lrgsgpath, 
         'ordParam': lattice.phtrapath}
filename = file_path_maker(mpath[mode])
if os.path.exists(filename):
    exit(f"File {os.path.split(filename)[1]} already exists.")
else:
    nAvgDone = 0
    try:
        fnameExists = glob.glob(f"{file_path_maker(mpath[mode], napath='', 
                                               spath='', extout='')}*")[0]
        merged_dict = pk.load(open(fnameExists, 'rb'))   
        if outsx:
            avgIdx = -2
        else:
            avgIdx = -1
        nAvgDone = os.path.splitext(fnameExists.split('_')[avgIdx])[0]
        nAvgDone = int(re.search(r'\d+', nAvgDone).group())
        fnameOld = fnameExists
    except:
        fnameOld = file_path_maker(mpath[mode], napath=0)
    nAvgNeed = navg - nAvgDone
#
#
#
match mode:
    case 'pCluster':
        period = sfreq
        for cp in range((nAvgNeed // period) + bool(nAvgNeed % period)):
            batch_size = min(nAvgNeed - cp * period, period)
            for _ in range(batch_size):
                l = Lattice3D(side, pflip=p, geo=geo, pdil=pdil, 
                    init_nw_dict=True)
                l.flip_sel_edges(geometry_func(l))
                #
                dist_dict = l.get_cluster_distribution()
                merged_dict += dist_dict
            navgCurr = nAvgDone + (cp + 1) * period
            fnameNew = file_path_maker(mpath[mode], napath=navgCurr)
            try:
                os.rename(fnameOld, fnameNew)
            except FileNotFoundError or OSError:
                pass
            with open(fnameNew, "wb") as f:
                pk.dump(merged_dict, f)
            fnameOld = fnameNew
    case 'ordParam':
        Pinf = []
        Pinf2 = []
        neglinks = 0
        for cont, avg in enumerate(range(navg)):
            avg1 = avg+1
            l = Lattice3D(side, pflip=p, geo=geo, pdil=pdil, 
                init_nw_dict=True)
            l.flip_sel_edges(geometry_func(l))
            neglinks += l.Ne_n
            #
            l.compute_k_eigvV(typf=typf)
            try:
                l.compute_pinf()
            except IndexError:
                continue
            #
            Pinf.append(l.Pinf)
            Pinf2.append(l.Pinf**2)
            #
            Pinf_tot = sum(Pinf)/avg1
            Pinf2_tot = sum(Pinf2)/avg1
            data=[avg1,
                l.pflip,
                neglinks/avg1,
                Pinf_tot,
                Pinf2_tot,
                Pinf2_tot-Pinf_tot**2]
            #
            if (avg1 % sfreq == 0):
                try:
                    filenameold = file_path_maker(mpath[mode], 
                                                  napath=avg1-sfreq)
                    os.remove(filenameold)
                except OSError:
                    pass
                filename = file_path_maker(mpath[mode], napath=avg1)
                with open(filename, 'wb') as file:
                    np.savetxt(file, np.atleast_2d(data), fmt='%.7g')