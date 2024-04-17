from parsers.Lattice2D_SlaplSpect_Serialiser_Parser import *

#

args = parser.parse_args()
#
geo = args.geometry
cell = args.cell_type
mode = args.mode
eigmode = args.eigen_mode
navg = args.number_of_averages
T = args.period
bins_count = args.bins_count
howmany = args.howmany
#
progName = Lattice2D_SlaplSpect_progName
progNameShrt = Lattice2D_SlaplSpect_progNameShrt
progMode = args.mode.split("_")[-1]
execBool = args.exec
printBool = args.print
#
if mode.endswith("eigDistr"):
    List = 2 ** np.arange(4, 10)
    plist = np.linspace(0.005, 0.5, num=10)
if mode.startswith("slanzarv"):
    if args.slanzarv_minMB == args.slanzarv_maxMB:

        def memoryfunc(*_):
            return args.slanzarv_minMB

    else:

        def memoryfunc(x):
            return int(
                np.interp(
                    x,
                    [min(List), max(List)],
                    [args.slanzarv_minMB, args.slanzarv_maxMB],
                )
            )

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
        argstr = f"{L} {p:.3g} -g {geo} -c {cell} -n {navg} --mode={mode}"
        return f"{slanzarv_str(mode, L, p, geo, cell)} {lnchStr} {argstr}"

else:
    exit(0)

count = 0
if mode.endswith("eigDistr"):
    for L in List:
        for p in plist:
            estring = exec_str(L, p, geo, cell, navg, progMode)
            operate(estring, count)