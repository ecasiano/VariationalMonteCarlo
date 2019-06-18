# Plots ground state energies of the BH model as a fucntion of variational parameters

import numpy as np
import matplotlib.pyplot as plt

#Load data
file_name = "Data/ddcerr_4_4_9.9763_1.6212_100000.dat"
file_name = "Data/ddcerr_10_10_9.9763_1.6212_100000.dat"
r = 1 #Set the distance between sites r
data = np.loadtxt(file_name)[:,r]  #The 
print(data)

#Extract BH parameters from file name
L,N,U,v,M = file_name.split("_")[1:]
M = M.split(".")[0]
L,N,U,v,M = int(L),int(N),float(U),float(v),int(M) #Promote from str to int OR float

#Determine the number of bin levels
bin_levels = [float(i) for i in range(np.size(data))]

#Plot
fig, ax1 = plt.subplots()
ax1.plot(bin_levels,data,'-',color='lightskyblue',label='9.9763')
ax1.plot(bin_levels,data,'o',color='steelblue',label='9.9763')
#ax1.axhline(y=-1.71173196633913,linewidth=1,color="#cccccc",zorder=0)
#ax1.text(1.5,-0.4,r'Exact Egs: $-1.71173196633913$')
ax1.set_ylabel(r"$Std. Error$")
ax1.set_xlabel(r"$Bin Level$")
#ax1.set_xlim(data00[:,0][0],data00[:,0][325])
#ax1.set_xlim(data16b[:,0][0],data16b[:,0][-1])
#ax1.set_ylim(-2,6)
#plt.legend(ncol=2,title=r"$U$")

plt.savefig("ddcerr_%i_%i_%.4f_%.4f_%i.pdf"%(L,N,U,v,M))
