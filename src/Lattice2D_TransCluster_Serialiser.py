import os
import sys
import numpy as np
#
programName = "Lattice2D_TransCluster"
programNamesht = "L2D_TC"
launchstr = f"python src/{programName}.py"
do_print = True
flag_mmemb = False
slanzarv_OPT = False
default_noAvg = 1000
mode = 'phaseTr'
List = 2**np.arange(4, 11)
plist = {L: np.concatenate(
            (
                np.linspace(1./L**1.5, 0.08, num=3),
                np.linspace(0.08, 0.2, num=10),
                np.linspace(0.2, 0.5, num=3)
            )
        ) for L in List}
if any(pmsg in sys.argv for pmsg in ["--pCluster", "-pC"]):
    mode = 'pCluster'
    plist = {L: np.concatenate(
            (
                np.linspace(1./L**1.5, 0.2, num=25),
                np.linspace(0.2, 0.5, num=10),
                np.linspace(0.5, 1, num=5)
            )
        ) for L in List}
if any(pmsg in sys.argv for pmsg in ["--execute", "-e"]):
    do_print = False
if any(pmsg in sys.argv for pmsg in ["--verbose", "-v", "--print"]):
    do_print = True
if "-n" in sys.argv:
    index = sys.argv.index("-n")
    noAvg = sys.argv[index+1]
else:
    noAvg = default_noAvg

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
    def slanzarv_STR(mode, L, p, geo, c):
        slanzarvstr = f'slanzarv -m {memoryfunc(L)} --nomail --jobname '
        argstr = f'"{programNamesht}{mode}_{L}_{p:.3g}_{geo}_{c}"'
        return slanzarvstr + argstr if slanzarv_OPT else ""        
else:
    def slanzarv_STR(*args):
        return ""
#
count = 0
if mode == 'pCluster':
    if "-g" in sys.argv:
        index = sys.argv.index("-g")
        geo = sys.argv[index+1]
    else:
        geo = 'squared'
    if "-c" in sys.argv:
        index = sys.argv.index("-g")
        cell = sys.argv[index+1]
    else:
        cell = 'single'
    #
    for L in List:
        for p in plist[L]:
            argstr = f"{L} {p:.3g} -g {geo} -c {cell} -nA {noAvg} --mode={mode}"
            the_string = f"{slanzarv_STR(mode, L, p, geo, cell)} {launchstr} {argstr}"
            if do_print:
                print(the_string)
            else:
                os.system(the_string)
            count += 1
elif mode == 'phaseTr':
    geometry_cell_dict = {'squared': ['single', 'singleZERR', 'singleXERR'],
                        'triangular': ['single', 'singleZERR', 'singleXERR'],
                        'hexagonal': ['single', 'singleZERR', 'singleXERR']}
    #
    for L in List:
        for p in plist[L]:
            for geo, cellst in geometry_cell_dict.items():
                for c in cellst:
                    argstr = f"{L} {p:.3g} -g {geo} -c {c} -nA {noAvg} --mode={mode}"
                    the_string = f"{slanzarv_STR(mode, L, p, geo, c)} {launchstr} {argstr}"
                    if do_print:
                        print(the_string)
                    else:
                        os.system(the_string)
                    count += 1
print("submitted jobs: ", count)