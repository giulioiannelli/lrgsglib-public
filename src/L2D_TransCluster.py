from parsers.L2D_TransCluster_Parser import *
#
args = parser.parse_args()
#
side = args.L
p = args.p
geo = args.geometry
cell = args.cell_type
mode = args.mode
navg = args.number_of_averages
sfreq = args.save_frequency if args.save_frequency else navg // 20
outsx = args.out_suffix
typf = args.float_type
#
match typf:
    case 'float32':
        typf = np.float32
    case 'float64':
        typf = np.float64
    case _:
        raise ValueError("Invalid float type specified")
#
def get_geometry_func(cell: str):
    match cell:
        case 'rand' | 'randZERR' | 'randXERR':
            def geometry_func(lattice: Lattice2D):
                return lattice.nwDict[cell]['G']
        case _ if cell.startswith('ball'):
            radius = get_first_int_in_str(cell)
            def geometry_func(lattice: Lattice2D):
                return lattice.nwDict.get_links_rball(radius)
        case _:
            raise ValueError("Invalid cell specified")

    return geometry_func
geometry_func = get_geometry_func(cell)
# if cell in ['rand', 'randZERR', 'randXERR'] or cell.startswith('ball'):
#     if cell == 'rand':
#         def geometry_func(lattice: Lattice2D):
#             return lattice.nwDict[cell]['G']
#     elif cell == 'randZERR':
#         def geometry_func(lattice: Lattice2D):
#             return lattice.nwDict[cell]['G']
#     elif cell == 'randXERR':
#         def geometry_func(lattice: Lattice2D):
#             return lattice.nwDict[cell]['G']
#     elif cell.startswith('ball'):
#         radius = get_first_int_in_str(cell)
#         def geometry_func(lattice: Lattice2D):
#             return lattice.nwDict.get_links_rball(radius)
# else:
#     raise ValueError("Invalid cell specified")
#
if mode == 'pCluster':
    merged_dict = Counter()
    extout = PKL
elif mode == 'ordParam':
    Pinf = np.zeros(navg)
    Pinf2 = np.zeros(navg)
    Fluct = np.zeros(navg)
    Fluct2 = np.zeros(navg)
    extout = TXT
else:
    raise ValueError("Invalid mode specified")
#
def file_path_maker(mpath, mode=mode, ppath = p, 
                    napath = navg, 
                    spath = outsx,
                    ctpath = cell,
                    extout = extout):
    if spath:
        spath = "_"+spath
    return f'{mpath}{mode}_p={ppath:.3g}_{ctpath}_na={napath}{spath}{extout}'
#
lattice = Lattice2D(side, pflip=p, geo=geo)
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
if mode == 'pCluster':
    period = sfreq
    for current_period in range((nAvgNeed // period) + bool(nAvgNeed % period)):
        batch_size = min(nAvgNeed - current_period * period, period)
        for _ in range(batch_size):
            l = Lattice2D(side, pflip=p, geo=geo, initNwDict=True)
            l.flip_sel_edges(geometry_func(l))
            #
            l.compute_k_eigvV(typf=typf)
            dist_dict = l.cluster_distribution()
            merged_dict += dist_dict
        navgCurr = nAvgDone + (current_period + 1) * period
        fnameNew = file_path_maker(mpath[mode], napath=navgCurr)
        try:
            os.rename(fnameOld, fnameNew)
        except FileNotFoundError or OSError:
            pass
        with open(fnameNew, "wb") as f:
            pk.dump(merged_dict, f)
        fnameOld = fnameNew
elif mode == 'ordParam':
    neglinks = 0
    for cont, avg in enumerate(range(navg)):
        avg1 = avg+1
        l = Lattice2D(side, pflip=p, geo=geo, initNwDict=True)
        l.flip_sel_edges(geometry_func(l))
        #
        l.compute_k_eigvV(typf=typf)
        l.calc_Pinf()
        #
        Pinf[cont]=l.Pinf
        Pinf2[cont]=l.Pinf**2
        #
        neglinks += l.Ne_n
        Pinf_tot = np.sum(Pinf)/avg1
        Pinf2_tot = np.sum(Pinf2)/avg1
        data=[avg1, 
              l.pflip,
              neglinks/avg1,
              Pinf_tot, 
              Pinf2_tot,
              Pinf_tot**2-Pinf2_tot]
        #
        if (avg % sfreq == 0):
            try:
                filenameold = file_path_maker(mpath[mode], napath=avg)
                os.remove(filenameold)
            except OSError:
                pass
            filename = file_path_maker(mpath[mode], napath=avg+sfreq)
            with open(filename, 'wb') as file:
                np.savetxt(file, np.atleast_2d(data), fmt='%.7g')