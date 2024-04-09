import numpy as np
import os

s = [5, 10, 20, 30, 40, 50]
p = np.linspace(0.15, 0.25, num=10)
printb = True

for side in s:
    for pflip in p:
        execstr = f"slanzarv -m 16384 --nomail python src/Lattice3D_TransCluster.py " + \
    f"{side} {pflip:.3g} -n 10000 --mode=pCluster"
        if printb:
            print(execstr)
        else:
            os.system(execstr)