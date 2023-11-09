import numpy as np
import os
import sys

L = int(sys.argv[1])
p = float(sys.argv[2])
T_l = np.linspace(0.005, 0.5, 100)

for T in T_l:
    os.system(
        "slanzarv python3.10 lattc2dsq_IsingDynamics_mClusters.py "
        + f"{L} {p:.3g} {T:.4g} -nA 1000"
    )
