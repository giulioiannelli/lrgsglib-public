from parsers.L2D_IsingDynamics_Serialiser import parse_arguments, parser
from lrgsglib.config.progargs import *
from lrgsglib.shared import *
#
args = parse_arguments(parser)
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
thrmsteps = args.thermsteps
ic_gs = args.init_cond.startswith('ground_state')
number = int(args.init_cond.split('_')[-1]) if ic_gs else 0
workdir = args.workdir
#
progName = L2D_IsingDynamics_progName
progNameShrt = L2D_IsingDynamics_progNameShrt
execBool = args.exec
printBool = args.print
#
List = args.side1_list
navglist = [navg for _ in range(len(List))]
plist = args.pflip_linsp
Tlist = args.Temp_linsp

if args.slanzarv_minMB == args.slanzarv_maxMB:
    def memoryfunc(*_):
        return args.slanzarv_minMB
else:
    def memoryfunc(x):
        return np.interp(x, [min(List), max(List)],
                [args.slanzarv_minMB, args.slanzarv_maxMB]).astype(int)
def slanzarv_str(mode, L, p, geo, c, T, nomail=True, moretime=False, short=True):
    jobname = '_'.join([progNameShrt, mode[:3], f"{L}", f"{p:.3g}", f"{T:.3g}", 
                        geo[:3], out_suffix])
    slanzarvopt = f"--jobname {jobname} "
    if nomail:
        slanzarvopt += "--nomail "
    if moretime:
        slanzarvopt += " --time " + moretime + ""
    if short:
        slanzarvopt += "--short"
    slanzarvstr = f"slanzarv -m {memoryfunc(L)} {slanzarvopt}"
    return slanzarvstr
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
                 NoClust, workdir=workdir):
        lnchStr = f"python src/{progName}.py"
        in_suffix = in_suffix or f"p={p:.3g}"
        outstr = f"-os {out_suffix} " if out_suffix else " "
        argstr = (f"{L} {p:.3g} {T:.3g} -g {geo} -c {cell} -n {navg} -ic {ic} "
                  f"-rl {runlang} -is {in_suffix} "
                  f"{outstr} "
                  f"-nc {NoClust} "
                  f"-ts {thrmsteps} "
                  f"-wd {workdir} " if workdir else "")
        return f"{slanzarv_str(runlang, L, p, geo, cell, T)} {lnchStr} {argstr}"

else:
    exit(0)

for L,navv in zip(List, navglist):
    for p in plist:
            for T in Tlist:
                estring = exec_str(L, p, geo, cell, navv, T, 
                                   runlang, in_suffix, out_suffix, NoClust)
                operate(estring)
print(f"Total number of jobs executed: {count_exe}")
print(f"Total number of jobs printed: {count}")