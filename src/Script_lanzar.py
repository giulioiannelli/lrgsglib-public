import numpy as np
import os
import sys 

L=int(sys.argv[1])
p2=float(sys.argv[2])
t1=np.linspace(0.005,0.5,100)


for t in t1:
	#print('slanzarv ./VortexProg 10 '+str(t)+' 1.0 AdyNet'+str(nl)+' AdyNetX'+str(nl)+' 310')
	os.system('slanzarv python3.10 lattc2dsq_IsingDynamics2.py '+str(L)+' 4 '+str(p2)+' '+str(t)+' -nA 1000')
