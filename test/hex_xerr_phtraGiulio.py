import subprocess
import numpy as np

script_path = "src/L2D_TransCluster.py"

for p in np.linspace(0, 0.5, 20):
	args = ["32", f"{p:.3g}", "-c", "randXERR", "-n", "50", "--outdir", 
		"data/test_hexerr/", "-g", "hex"]
	subprocess.run(["python", script_path] + args)
