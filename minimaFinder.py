#Do a quadratic fitting of given array of data
#then find the minimum value of the fit (vertex)
import numpy as np

data = np.loadtxt('Data/EOPP6F3l3a2.dat')

quadFit = np.polyfit(data[:,0],data[:,3],deg=2)

uMax = -(quadFit[1]/(2*quadFit[0]))
sMax = quadFit[0]*uMax**2 + quadFit[1]*uMax + quadFit[2]

print('Max: V/T = %.16f, S1_OP = %.16f'%(uMax,sMax))