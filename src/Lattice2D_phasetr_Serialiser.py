import numpy as np
import os
import sys
#
List = [256, 512]
plist = np.concatenate((np.logspace(-3, np.log10(0.05), num=5),
                        np.linspace(0.05, 0.2, num=10),
                        np.linspace(0.2, 0.5, num=5)))
geometry_cell_dict = {'squared': ['single', 'square', 'cross'],
                      'triangular': ['single', 'triangle', 'cross'],
                      'hexagonal': ['single', 'hexagon', 'cross']}

programName = "Lattice2D_phasetr"
launchstr = f"python src/{programName}.py"
do_print = True
flag_mmemb = False
slanzarv_OPT = False
default_noAvg = 1000
#
if any(pmsg in sys.argv for pmsg in ["--execute", "-e"]):
    do_print = False
if any(pmsg in sys.argv for pmsg in ["--verbose", "-v", "--print"]):
    do_print = True
#
if "-n" in sys.argv:
    index = sys.argv.index("-n")
    noAvg = sys.argv[index+1]
else:
    noAvg = default_noAvg
#
if "--slanzarv" in sys.argv:
    slanzarv_OPT = True
    index = sys.argv.index("--slanzarv")
    try:
        slanzarv_mEMMB = int(sys.argv[index+1])
    except IndexError:
        slanzarv_mEMMB = 1024
    try:
        slanzarv_MEMMB = int(sys.argv[index+2])
    except IndexError:
        slanzarv_MEMMB = 10240
        flag_mmemb =True
    if flag_mmemb:
        def memoryfunc(x):
            return slanzarv_mEMMB
    else:
        def linear_map(x, x1=min(List), x2=max(List), y1=slanzarv_mEMMB, 
                    y2=slanzarv_MEMMB):
            a = (y2 - y1) / (x2 - x1)
            b = y1 - a * x1
            return a * x + b
        def memoryfunc(x):
            return int(linear_map(x))
    def slanzarv_STR(L, p):
        slanzarvstr = f'slanzarv -m {memoryfunc(L)} --nomail --jobname'
        argstr = f'"{programName}_{L}_{p:.3g}"'
        return slanzarvstr + argstr if slanzarv_OPT else ""
            
else:
    def slanzarv_STR(*args):
        return ""
#
#
#
for L in List:
    for p in plist:
        for geo, cellst in geometry_cell_dict.items():
            for c in cellst:
                argstr = f"{L} {p:.3g} -g {geo} -c {c} -nA {noAvg}"
                the_string = f"{slanzarv_STR(L, p)} {launchstr} {argstr}"
                print(the_string)
                if not do_print:
                    os.system(the_string)
