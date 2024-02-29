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
do_print = sys.argv[1]
slanzarv_OPT = sys.argv[2]
slanzarv_MEMMB = sys.argv[3]

progname = "lattc2d_phasetr"

for L in List:
    for p in plist:
        for geo, cellst in geometry_cell_dict.items():
            for c in cellst:
                slanzarv_STR = (
                    f'slanzarv -m {slanzarv_MEMMB} --nomail --jobname "{progname}_{L}_{p:.3g}"'
                    if slanzarv_OPT == "--slanzarv"
                    else ""
                )
                the_string = (
                    slanzarv_STR
                    + f" python src/{progname}.py "
                    + f"{L} {p:.3g} -g {geo} -c {c} -nA {noavg} "
                )
                if do_print == "--print":
                    print(the_string)
                else:
                    print(the_string)
                    os.system(the_string)
