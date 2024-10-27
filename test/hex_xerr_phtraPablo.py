import subprocess
import numpy as np

script_path = "src/Code_PhT_Glass.py"

for p in np.linspace(0, 0.5, 20):
	args = ["32", f"{p:.3g}"]
	subprocess.run(["python", script_path] + args)
