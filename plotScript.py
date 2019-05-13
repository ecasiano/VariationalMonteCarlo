# Plots ground state energies of the BH model as a fucntion of variational parameters

import numpy as np
import matplotlib.pyplot as plt

#Load data
data00 = np.loadtxt("Data/EGS_4_4_0.0000.dat")
data01 = np.loadtxt("Data/EGS_4_4_0.0500.dat")
data02 = np.loadtxt("Data/EGS_4_4_0.2500.dat")
data03 = np.loadtxt("Data/EGS_4_4_0.5000.dat")
data04 = np.loadtxt("Data/EGS_4_4_0.7500.dat")
data05 = np.loadtxt("Data/EGS_4_4_1.0000.dat")
data06 = np.loadtxt("Data/EGS_4_4_1.2500.dat")
data07 = np.loadtxt("Data/EGS_4_4_1.5000.dat")
data08 = np.loadtxt("Data/EGS_4_4_1.7500.dat")
data09 = np.loadtxt("Data/EGS_4_4_2.0000.dat")
data10 = np.loadtxt("Data/EGS_4_4_2.2500.dat")
data11 = np.loadtxt("Data/EGS_4_4_2.5000.dat")
data12 = np.loadtxt("Data/EGS_4_4_2.7500.dat")
data13 = np.loadtxt("Data/EGS_4_4_3.0000.dat")
data14 = np.loadtxt("Data/EGS_4_4_3.2500.dat")
data15 = np.loadtxt("Data/EGS_4_4_3.5000.dat")


#data05 = np.loadtxt("EGS_4_4_199.0536.dat")

#Plot
fig, ax1 = plt.subplots()
ax1.plot(data15[:,0],data15[:,1],'-',label='3.5000')
ax1.plot(data14[:,0],data14[:,1],'-',label='3.2500')
ax1.plot(data13[:,0],data13[:,1],'-',label='3.0000')
ax1.plot(data12[:,0],data12[:,1],'-',label='2.7500')
ax1.plot(data11[:,0],data11[:,1],'-',label='2.5000')
ax1.plot(data10[:,0],data10[:,1],'-',label='2.2500')
ax1.plot(data09[:,0],data09[:,1],'-',label='2.0000')
ax1.plot(data08[:,0],data08[:,1],'-',label='1.7500')
ax1.plot(data07[:,0],data07[:,1],'-',label='1.5000')
ax1.plot(data06[:,0],data06[:,1],'-',label='1.2500')
ax1.plot(data05[:,0],data05[:,1],'-',label='1.0000')
ax1.plot(data04[:,0],data04[:,1],'-',label='0.7500')
ax1.plot(data03[:,0],data03[:,1],'-',label='0.5000')
ax1.plot(data02[:,0],data02[:,1],'-',label='0.2500')
ax1.plot(data01[:,0],data01[:,1],'-',label='0.0500')
ax1.plot(data00[:,0],data00[:,1],'-',label='0.0000')
#ax1.plot(data04[:,0],data04[:,1],'-',label='v=5.1165')
#ax1.plot(data05[:,0],data05[:,1],'-',label='v=7.0627')
#ax1.plot(data06[:,0],data06[:,1],'-',label='v=31.5479')
#ax1.plot(data02[:,0],data02[:,1],'-',label='v=5.1165')
#ax1.plot(data03[:,0],data03[:,1],'-',label='v=199.0536')
ax1.set_ylabel(r"$Egs$")
ax1.set_xlabel(r"$v$")
ax1.set_xlim(data00[:,0][0],data00[:,0][325])
plt.legend(ncol=2,title=r"$U$")
plt.savefig('groundStatesVariousParams.pdf')