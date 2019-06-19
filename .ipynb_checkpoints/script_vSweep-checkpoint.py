#
import subprocess
import numpy as np

program = 'boseHubbardVMC.py'
L = '4'
N = '4'
U = '9.97631157484440'
v_min = 1.25
v_max = 2.001
v_step = 0.002
v_list = np.arange(v_min,v_max,v_step)
v_list = np.arange(1.25,2.001,0.01)
M = str(int(1.5E+05))
t = '1.0'
mu = '0.0'

for v in v_list:
    v = str(v)
    subprocess.call(['python', program, '--t', t, '--mu', mu, L, N, U, v, M])                