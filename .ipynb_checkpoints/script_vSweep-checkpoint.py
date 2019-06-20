#
import subprocess
import numpy as np

program = 'boseHubbardVMC.py'
L = '4'
N = '4'
U = '9.97631157484440'
v_min = 1.36
v_max = 2.0001
v_step = 0.001
v_list = np.arange(v_min,v_max,v_step)
M = str(int(1.5E+05))
t = '1.0'
mu = '0.0'

print("Completion: ")
for i,v in enumerate(v_list):
    v = str(v)
    subprocess.call(['python', program, '--t', t, '--mu', mu, L, N, U, v, M])  
    print("%.2f %%"%((i+1)/np.shape(v_list)[0]*100))
