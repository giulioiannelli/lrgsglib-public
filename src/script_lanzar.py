import numpy as np
import os
import sys

L = int(sys.argv[1])
p2 = float(sys.argv[2])
t1 = np.linspace(0.005, 0.5, 100)

for t in t1:
    os.system(
        "slanzarv python3.10 lattc2dsq_IsingDynamics_mClusters.py "
        + str(L)
        + " "
        + str(p2)
        + " "
        + str(t)
        + " -nA 1000"
    )
