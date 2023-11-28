import numpy as np
import os
import sys

L = int(sys.argv[1])
p_l = np.linspace(0.075, 0.15, 10)#float(sys.argv[2])
T_l = np.linspace(0.005, 0.5, 100)

for p in p_l:
    for T in T_l:
        os.system(
            "slanzarv python lattc2dsq_IsingDynamics_mClusters.py "
            + f"{L} {p:.3g} {T:.4g} -nA 1000"
        )
        # print(
        #     "slanzarv python lattc2dsq_IsingDynamics_mClusters.py "
        #     + f"{L} {p:.3g} {T:.4g} -nA 1000"
        # )

