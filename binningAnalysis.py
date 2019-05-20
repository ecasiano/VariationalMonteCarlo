#Perform binning analysis of observables obtained
#via boseHubbardVMC.py

import numpy as np
import matplotlib.pyplot as plt

def get_binned_error(mc_data):
    '''Get the standard error in mc_data and return neighbor averaged data.'''
    N_bins = mc_data.size
    std_error = np.std(mc_data)/np.sqrt(N_bins)     #Standard error
    
    start_bin = N_bins % 2
    binned_mc_data = 0.5*(mc_data[start_bin::2]+mc_data[start_bin+1::2]) #Averages (A0,A1), (A2,A3), + ... A0 ignored if odd data
   
    return std_error,binned_mc_data

def main():
    
    #Load data
    data = np.loadtxt("Data/EGS_4_4_9.9763.dat")
    
    
    ##### NEED: 1. Find the optimal variational parameter 2. Use it to generate egs vs mc_step data for binning analysis
    
    
    #delta, binned_egs = get_binned_error(data[])
    
if __name__ == "__main__":
    main()