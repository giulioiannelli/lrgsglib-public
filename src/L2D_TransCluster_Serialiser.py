from parsers.L2D_TransCluster_Serialiser_Parser import *
#
args = parser.parse_args()
#
progName = Lattice2D_TransCluster_progName
progNameShrt = Lattice2D_TransCluster_progNameShrt
navg = args.number_of_averages
fullMode = args.mode
progMode = args.mode.split('_')[-1]
execBool = args.exec
printBool = args.print
#
if fullMode.endswith('pCluster'):
    List = [16, 32, 64, 96, 128, 192, 256, 512]
    geo = args.geometry
    cell = args.cell_type
    plist = np.linspace(0.09, .2, num=20)
elif fullMode.endswith('ordParam'):
    List = [192, 256, 512]
    def linspacepfunc(geo, cell):
        if cell != 'randXERR':
            if geo != 'hexagonal':
                return np.linspace(0, 0.3, 100)
            else:
                return np.linspace(0, 0.35, 100)
        return np.linspace(0, 1, 100)
    geometry_cell_dict = {'squared': L2D_RAND_CELL_LIST,
                        'triangular': L2D_RAND_CELL_LIST,
                        'hexagonal': L2D_RAND_CELL_LIST}
    plist = {geo: {cell: linspacepfunc(geo, cell) for cell in cells}  
            for geo,cells in geometry_cell_dict.items()}
#
if fullMode.startswith('slanzarv'):
    if args.slanzarv_minMB == args.slanzarv_maxMB:
        def memoryfunc(*_):
            return args.slanzarv_minMB
    else:
        def memoryfunc(x):
            return int(np.interp(x, [min(List), max(List)], 
                             [args.slanzarv_minMB, args.slanzarv_maxMB]))
    def slanzarv_str(mode, L, p, geo, c):
        slanzarvopt = "--nomail --jobname "
        slanzarvstr = f"slanzarv -m {memoryfunc(L)} {slanzarvopt}"
        argstr = f"{progNameShrt}{mode}_{L}_{p:.3g}_{geo[:3]}_{c[3:]}"
        return slanzarvstr + argstr 
else:
    def slanzarv_str(*_):
        return ""
#
if execBool or printBool:
    if execBool and printBool:
        def operate(s, count):
            print(s)
            os.system(s)
            count += 1
    elif execBool:
        def operate(s, count):
            os.system(s)
            count += 1
    elif printBool:
        def operate(s, *args):
            print(s)
    def exec_str(L, p, geo, cell, navg, mode):
        lnchStr = f"python src/{progName}.py"
        argstr = (f"{L} {p:.3g} -g {geo} -c {cell} -n "
                    f"{navg} --mode={mode}")
        return (f"{slanzarv_str(mode, L, p, geo, cell)} "
                        f"{lnchStr} {argstr}")
else:
    exit(0)

count = 0
if fullMode.endswith('pCluster'):
    for L in List:
        for p in plist:
            estring = exec_str(L, p, geo, cell, navg, progMode)
            operate(estring, count)
elif fullMode.endswith('ordParam'):
    for L in List:
        for geo, cellst in geometry_cell_dict.items():
            for cell in cellst:
                for p in plist[geo][cell]:
                    estring = exec_str(L, p, geo, cell, navg, progMode)
                    operate(estring, count)
else:
    print(f"no program executed, unkonw mode provided: {fullMode}")
    exit(0)         
print("submitted jobs: ", count)