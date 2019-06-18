#Perform binning analysis of observables obtained
#via boseHubbardVMC.py

import numpy as np
import matplotlib.pyplot as plt

def get_std_error(mc_data):
    '''Input array and calculate standard error'''
    N_bins = np.shape(mc_data)[0]
    std_error = np.std(mc_data)/np.sqrt(N_bins)
    
    return std_error

def get_binned_data(mc_data):
    '''Return neighbor averaged data.'''
    N_bins = np.shape(mc_data)[0]
    start_bin = N_bins % 2
    binned_mc_data = 0.5*(mc_data[start_bin::2]+mc_data[start_bin+1::2]) #Averages (A0,A1), (A2,A3), + ... A0 ignored if odd data

    return binned_mc_data

def get_autocorrelation_time(error_data):
    '''Given an array of standard errors, calculates autocorrelation time'''
    print(error_data[0],error_data[-2])
    autocorr_time = 0.5*((error_data[-2]/error_data[0])**2 - 1)
    return autocorr_time

'''------------------------------------------------------------------------'''

def main():

    #Load data
    file_name = "Data/ddc_4_4_9.9763_1.6212_100000.dat" 
    file_name = "Data/ddc_10_10_9.9763_1.6212_100000.dat"
    data = np.loadtxt(file_name)
    
    #Set from where to start the data depending on equilibration time
    equil_time = 0.2 #Set the percentage of unequilibrated data to throw away
    begin_data = int(np.shape(data)[0]*equil_time)
    data = np.loadtxt(file_name)[begin_data:]

    #Extract BH parameters from file name
    L,N,U,v,M = file_name.split("_")[1:]
    M = M.split(".")[0]
    L,N,U,v,M = int(L),int(N),float(U),float(v),int(M) #Promote from str to int OR float
    
    #Determine max bin level
    max_bin_level = int(np.log2(np.shape(data)[0]))
    min_bin = 40
    
    #Initialize list to save standard error
    std_errors = []
    r_max = np.shape(data)[1]
    
    #Binning loop
    binned_data = np.copy(data)
    for i in range(max_bin_level):
        std_errors_i = []
        print(np.shape(binned_data)[0])
        for r in range(r_max):
            std_errors_i.append(get_std_error(binned_data[:,r])) #error for each length r at bin level i
        std_errors.append(std_errors_i)   
        if np.shape(binned_data)[0]/2 <= min_bin: break
        binned_data = get_binned_data(binned_data)
                            
    #print("<E_gs>: (Raw) %.12f (Binned) %.12f"%(np.mean(original_data),np.mean(binned_data)))
    
    #Format the data file
    std_errors = np.array(std_errors)
    with open("ddcerr_%i_%i_%.4f_%.4f_%i.dat"%(L,N,U,v,M),"w+") as data:
        np.savetxt(data,std_errors,delimiter=" ",fmt="%.16f",header="rows: mc_step columns: distance between sites (r=0,1,2,...) // \nBH Parameters: L=%d,N=%d,U=%.4f,v=%.14f,MC_steps=%d (%d after equilibration) "%(L,N,U,v,M,M*(1-equil_time)))
    

if __name__ == "__main__":
    main()
