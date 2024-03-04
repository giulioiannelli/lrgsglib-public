import numpy as np
import os
import sys
#
List = [16, 32]#64, 128, 256, 512, 1024
plist = np.concatenate((np.logspace(-3, np.log10(0.05), num=5),
                        np.linspace(0.05, 0.2, num=15),
                        np.linspace(0.2, 0.5, num=5)))
geometry_cell_dict = {'squared': ['single', 'square', 'cross'],
                      'triangular': ['single', 'triangle', 'cross'],
                      'hexagonal': ['single', 'hexagon', 'cross']}
noavg = 1000

programName = "Lattice2D_phasetr"

do_print = False
if "--print" in sys.argv:
    do_print = True
slanzarv_OPT = False
if "--slanzarv" in sys.argv:
    slanzarv_OPT = True
    try:
        index = sys.argv.index("--slanzarv")
        slanzarv_MEMMB = sys.argv[index+1]
    except IndexError:
        slanzarv_MEMMB = 1024
if "-n" in sys.argv:
    index = sys.argv.index("-n")
    noAvg = sys.argv[index+1]
else:
    noAvg = 1000
if "-g" in sys.argv:
    index = sys.argv.index("-g")
    geometry = sys.argv[index+1]
else:
    geometry = 'squared'

for L in List:
    for p in plist:
        for geo, cellst in geometry_cell_dict.items():
            for c in cellst:
                slanzarv_STR = (
                    f'slanzarv -m {slanzarv_MEMMB} --nomail --jobname "{programName}_{L}_{p:.3g}"'
                    if slanzarv_OPT == "--slanzarv"
                    else ""
                )
                the_string = (
                    slanzarv_STR
                    + f" python src/{programName}.py "
                    + f"{L} {p:.3g} -g {geo} -c {c} -nA {noavg} "
                )
                if do_print:
                    print(the_string)
                else:
                    print(the_string)
                    os.system(the_string)
