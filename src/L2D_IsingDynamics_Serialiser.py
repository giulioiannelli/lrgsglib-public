from parsers.L2D_IsingDynamics_Serialiser import parse_arguments, parser
from lrgsglib.config.funcs import *
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
progn = L2D_IsingDynamics_progName
progn_shrt = L2D_IsingDynamics_progNameShrt
exec_bool = args.exec
prnt_bool = args.print
#
side_list = args.side1_list
navg_list = [navg for _ in range(len(side_list))]
pflp_list = args.pflip_linsp
temp_list = args.Temp_linsp

if args.slanzarv_minMB == args.slanzarv_maxMB:
    def memoryfunc(*_):
        return args.slanzarv_minMB
else:
    def memoryfunc(x):
        hl_side = [min(side_list), max(side_list)]
        hl_memy = [args.slanzarv_minMB, args.slanzarv_maxMB]
        return np.interp(x, hl_side, hl_memy).astype(int)
def slanzarv_str(mode, L, p, geo, c, T, 
                 nomail=True, moretime=False, short=True):
    jobname = join_non_empty('_', progn_shrt, mode[:3],
                            f"{L}", f"{p:.3g}", f"{T:.3g}", geo[:3], c, 
                            out_suffix)
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
if exec_bool or prnt_bool:
    count = 0
    count_exe = 0
    def operate(s):
        global count, count_exe
        if exec_bool and prnt_bool:
            print(s)
            os.system(s)
            count += 1
            count_exe += 1
        elif exec_bool:
            os.system(s)
            count_exe += 1
        elif prnt_bool:
            print(s)
            count += 1
    def exec_str(L, p, geo, cell, navg, T, runlang, in_suffix, out_suffix, ic,
                 NoClust, workdir=workdir):
        lnchStr = f"python src/{progn}.py"
        instr = ' '.join([f"-is ", f"{in_suffix}" or f"p={p:.3g}"])
        outstr = f"-os {out_suffix} " if out_suffix else ""
        wdstr = f"-wd {workdir} " if workdir else ""
        argstr = join_non_empty(' ', f"{L}", f"{p:.3g}", f"{T:.3g}", 
                                f"-g {geo}", f"-c {cell}", f"-n {navg}", 
                                f"-ic {ic}", f"-rl {runlang}", instr,
                                outstr, f"-nc {NoClust}", f"-ts {thrmsteps}", wdstr)
        slanzarvStr = slanzarv_str(runlang, L, p, geo, cell, T)
        print(lnchStr, argstr, slanzarvStr)
        return ' '.join([slanzarvStr, lnchStr, argstr])

else:
    exit(0)

for L,navv in zip(side_list, navg_list):
    for p in pflp_list:
            for T in temp_list:
                estring = exec_str(L, p, geo, cell, navv, T, 
                                   runlang, in_suffix, out_suffix, ic, NoClust)
                operate(estring)
print(f"Total number of jobs executed: {count_exe}")
print(f"Total number of jobs printed: {count}")