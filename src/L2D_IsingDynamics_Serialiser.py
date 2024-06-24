from parsers.L2D_IsingDynamics_Serialiser import *
#
args = parser.parse_args()
#
geo = args.geometry
cell = args.cell_type
navg = args.number_of_averages
ic = args.init_cond
navg = args.number_of_averages
out_suffix = args.out_suffix
in_suffix = args.in_suffix
NoClust = args.NoClust
runlang = args.runlang
if ic.startswith('ground_state'):
    parts = ic.split('_')
    number = int(parts[-1])
    howmany = number+1
else:
    howmany = 1
#
progName = L2D_IsingDynamics_progName
progNameShrt = L2D_IsingDynamics_progNameShrt
execBool = args.exec
printBool = args.print
#
List = [16, 32, 64, 96]
plist = np.linspace(0.09, 0.4, num=10)
Tlist = np.concatenate([np.linspace(0.1, 1.4, 14),
    np.linspace(1.4, 2.5, 11),
    np.linspace(2.5, 5, 3)])

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
def slanzarv_str(mode, L, p, geo, c, T):
    slanzarvopt = "--nomail --jobname "
    slanzarvstr = f"slanzarv -m {memoryfunc(L)} {slanzarvopt}"
    argstr = '_'.join(progNameShrt, mode, f"{L}", f"{p:.3g}", f"{T:.3g}", 
                      geo[:3], c[3:])
    return slanzarvstr + argstr
#
if execBool or printBool:
    count = 0
    count_exe = 0
    def operate(s):
        global count, count_exe
        if execBool and printBool:
            print(s)
            os.system(s)
            count += 1
            count_exe += 1
        elif execBool:
            os.system(s)
            count_exe += 1
        elif printBool:
            print(s)
            count += 1
    def exec_str(L, p, geo, cell, navg, T, runlang, in_suffix, out_suffix, 
                 NoClust):
        lnchStr = f"python src/{progName}.py"
        in_suffix = in_suffix or f"p={p:.3g}"
        outstr = f"-os {out_suffix} " if out_suffix else " "
        argstr = (f"{L} {p:.3g} {T:.3g} -g {geo} -c {cell} -n {navg} -ic {ic} "
                  f"-rl {runlang} -is {in_suffix} "
                  f"{outstr}"
                  f"-nc {NoClust}")
        return f"{slanzarv_str(runlang, L, p, geo, cell, T)} {lnchStr} {argstr}"

else:
    exit(0)

for L in List:
    for p in plist:
            for T in Tlist:
                estring = exec_str(L, p, geo, cell, navg, T, 
                                   runlang, in_suffix, out_suffix, NoClust)
                operate(estring)
print(f"Total number of jobs executed: {count_exe}")
print(f"Total number of jobs printed: {count}")