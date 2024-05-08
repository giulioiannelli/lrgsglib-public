from parsers.L3D_TransCluster import *
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
#
if cell in ['rand', 'randZERR', 'randXERR'] or cell.startswith('ball'):
    if cell == 'rand':
        def geometry_func(lattice: Lattice3D):
            return lattice.nwDict[cell]['G']
    elif cell == 'randZERR':
        def geometry_func(lattice: Lattice3D):
            return lattice.nwDict[cell]['G']
    elif cell == 'randXERR':
        def geometry_func(lattice: Lattice3D):
            return lattice.nwDict[cell]['G']
    elif cell.startswith('ball'):
        radius = get_first_int_in_str(cell)
        def geometry_func(lattice: Lattice3D):
            return lattice.nwDict.get_links_rball(radius)
else:
    raise ValueError("Invalid cell specified")
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
lattice = Lattice3D(dim=(side for _ in range(3)), pflip=p, geo=geo)
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
            l = Lattice3D(side, pflip=p, geo=geo, initNwDict=True)
            l.flip_sel_edges(geometry_func(l))
            #
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
        l = Lattice3D(side, pflip=p, geo=geo, initNwDict=True)
        l.flip_sel_edges(geometry_func(l))
        #
        l.calc_fluct_Pinf()
        #
        Fluct[cont]=l.eigV_fluct
        Fluct2[cont]=l.eigV_fluct**2
        Pinf[cont]=l.Pinf
        Pinf2[cont]=l.Pinf**2
        #
        neglinks += l.Ne_n
        data=[l.pflip,
              neglinks,
              avg1,
              np.sum(Pinf)/avg1, 
              np.sum(Pinf2)/avg1, 
              np.sum(Fluct)/avg1, 
              np.sum(Fluct2)/avg1, 
              np.var(Fluct[Fluct!=0])]
        #
        if (avg % sfreq == 0):
            try:
                filenameold = file_path_maker(mpath[mode], napath=avg)
                os.remove(filenameold)
            except OSError:
                pass
            filename = file_path_maker(mpath[mode], napath=avg+sfreq)
            with open(filename, 'wb') as file:
                np.savetxt(file, data, fmt='%.7g')