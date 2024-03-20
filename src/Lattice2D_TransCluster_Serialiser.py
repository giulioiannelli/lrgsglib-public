from LRGSG_package.LRGSG import *
from parsers.Lattice2D_TransCluster_Serialiser_Parser import *
#
args = parser.parse_args()

progName = Lattice2D_TransCluster_progName
progNameShrt = f"{Lattice2D_TransCluster_progNameShrt}"
lnchStr = f"python src/{progName}.py"
navg = args.number_of_averages
# List = 2**np.arange(4, 9)
List = [8, 16, 32, 48, 64, 96, 128]
#
if args.mode.startswith('slanzarv'):
    if args.slanzarv_minMB == args.slanzarv_maxMB:
        def memoryfunc(x):
            return args.slanzarv_minMB
    else:
        def memoryfunc(x):
            return int(np.interp(x, [min(List), max(List)], 
                             [args.slanzarv_minMB, args.slanzarv_maxMB]))
                
    def slanzarv_STR(mode, L, p, geo, c):
        slanzarvstr = f'slanzarv -m {memoryfunc(L)} --nomail --jobname '
        argstr = f'"{progNameShrt}{mode}_{L}_{p:.3g}_{geo[:3]}_{c[4:]}"'
        return slanzarvstr + argstr 
else:
    def slanzarv_STR(*args):
        return ""
#
if args.exec or args.print:
    if args.exec and args.print:
        def operate(s, count):
            print(s)
            os.system(s)
            count += 1
    elif args.exec:
        def operate(s, count):
            os.system(s)
            count += 1
    elif args.print:
        def operate(s, *args):
            print(s)
    def exec_string(L, p, geo, cell, navg, mode):
        argstr = (f"{L} {p:.3g} -g {geo} -c {cell} -n "
                    f"{navg} --mode={mode}")
        return (f"{slanzarv_STR(mode, L, p, geo, cell)} "
                        f"{lnchStr} {argstr}")
    count = 0
    if args.mode.endswith('pCluster'):
        mode = 'pCluster'
        plist = {L: np.concatenate(
                (
                    np.linspace(1./L**1.5, 0.2, num=25),
                    np.linspace(0.2, 0.5, num=10),
                    np.linspace(0.5, 1, num=5)
                )
            ) for L in List}
        geo = args.geo
        cell = args.cell
        #
        for L in List:
            for p in plist[L]:
                argstr = f"""{L} {p:.3g} -g {geo} -c {cell} -A 
                    {args.number_of_averages} --mode={mode}"""
                the_string = f"""{slanzarv_STR(mode, L, p, geo, cell)} 
                    {lnchStr} {argstr}"""
                operate(the_string, count)
                count += 1
    elif args.mode.endswith('ordParam'):
        mode = 'ordParam'

        DEFLattice2D_singcellist = ['rand', 'randXERR', 'randZERR']
        geometry_cell_dict = {'squared': DEFLattice2D_singcellist,
                            'triangular': DEFLattice2D_singcellist,
                            'hexagonal': DEFLattice2D_singcellist}
        def linspacepfunc(g, c):
            if c != 'randXERR':
                if g != 'hexagonal':
                    return np.linspace(0, 0.3, 100)
                else:
                    return np.linspace(0, 0.35, 100)
            return np.linspace(0, 1, 100)

        plist = {geo: {cell: linspacepfunc(geo, cell) for cell in cells}  
                 for geo,cells in geometry_cell_dict.items()}
        #
        for L in List:
            for geo, cellst in geometry_cell_dict.items():
                for cell in cellst:
                    for p in plist[geo][cell]:
                        estring = exec_string(L, p, geo, cell, navg, mode)
                        operate(estring, count)         
    print("submitted jobs: ", count)