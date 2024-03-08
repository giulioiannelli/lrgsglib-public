import os
import sys
import numpy as np
#
programName = "Lattice2D_percClusters"
plist = [0.01, 0.034, 0.08, 0.103, 0.12, 0.206, 0.412]

List = [16, 32, 64, 128, 256, 512, 1024]
plist = {L: np.concatenate(
            (
                np.linspace(0.01, 0.09, num=3),
                np.linspace(0.09, 0.2, num=6),
                np.linspace(0.01, 0.09, num=3)
            )
        ) for L in List}


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
        slanzarv_STR = (
            f'slanzarv -m {slanzarv_MEMMB} --nomail --jobname "{programName}_{L}_{p:.3g}_{geometry}_{noAvg}"'
            if slanzarv_OPT
            else ''
        )
        the_string = (
            slanzarv_STR
            + f" python src/{programName}.py "
            + f"{L} {p:.3g} -g {geometry} -nA {noAvg} "
        )
        if do_print:
            print(the_string)
        else:
            print(the_string)
            os.system(the_string)
