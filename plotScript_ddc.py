# Plots ground state energies of the BH model as a fucntion of variational parameters

import numpy as np
import matplotlib.pyplot as plt

#Load data
file_name = "Data/ddc_10_10_9.9763_1.6212_100000.dat"
data = np.loadtxt(file_name)
errordata = np.loadtxt("Data/ddcerr_10_10_9.9763_1.6212_100000.dat")

#Set from where to start the data depending on equilibration time
equil_time = 0.2 #Set the percentage of unequilibrated data to throw away
begin_data = int(np.shape(data)[0]*equil_time)
data = np.loadtxt(file_name)[begin_data:,:]
#data = np.ones(10000)*2

#Extract BH parameters from file name
L,N,U,v,M = file_name.split("_")[1:]
M = M.split(".")[0]
L,N,U,v,M = int(L),int(N),float(U),float(v),int(M) #Promote from str to int OR float

#Calculate the average correlation function from raw data
r_max = np.size(data[0])
r_list = np.zeros(r_max)   #Store possible separations
c_list = np.zeros(r_max)   #Store average correlation functions
for r in range(r_max):
    c_list[r] = np.mean(data[:,r])
    r_list[r] = r
    
#Extract errorbars from ddc error files
r_max = np.shape(errordata)[1]
c_list_errors = np.zeros(r_max)
for r in range(r_max):
    c_list_errors[r] = np.max(errordata[:,r])

#Plot
fig, ax1 = plt.subplots()
ax1.errorbar(r_list, c_list, marker='o',ms=2.25,yerr=c_list_errors,fmt='o',capthick=10,
             ecolor='orange')
ax1.plot(r_list,c_list,'-',color='lightskyblue',label='9.9763',linewidth=0.25)
ax1.set_xticks(range(r_max))
ax1.set_ylabel(r"$\sum_i \langle n_i n_{i+r}\rangle$/L")
ax1.set_xlabel(r"$r$")
#ax1.set_xlim(data00[:,0][0],data00[:,0][325])
#ax1.set_xlim(data16b[:,0][0],data16b[:,0][-1])
#ax1.set_ylim(-2,6)
#plt.legend(ncol=2,title=r"$U$")
plt.savefig("ddc_%i_%i_%.4f_%.4f_%i.pdf"%(L,N,U,v,M))
