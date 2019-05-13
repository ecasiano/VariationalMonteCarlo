# Estimates the ground state energy for a Bose-Hubbard lattice
# using Variational Monte Carlo with a Jastrow Trial Wavefunction

import numpy as np
from sympy.utilities.iterables import multiset_permutations,ordered_partitions
import matplotlib.pyplot as plt

# Main function

# Function definitions

def bosonic_configurations(L,N):
    '''Input: 1D Lattice Size and Number of Bosons
    Output: All possible configurations of bosons'''

    #List that will store all configurations
    configurations = []

    #Store ordered partitions of N as a list
    partitions = list(ordered_partitions(N))

    for p in partitions:
        #BH Lattice containing a partition of N followed by zeros
        auxConfig = [0]*L
        auxConfig[0:len(p)] = p

        #Generate permutations based on current partition of N
        partitionConfigs = list(multiset_permutations(auxConfig))

        #Append permutations of current partition to list containing all configurations
        configurations += partitionConfigs

    #Promote configurations list to numpy array
    configurations = np.array(configurations)

    return configurations

'''------------------------------------------------------------------------'''

def random_boson_config(L,N):
    '''Generates a random configuration of N bosons in a 1D lattice of size L'''

    psi = np.zeros(L,dtype=int) # Stores the random configuration of bosons
    for i in range(N):
        r = np.random.randint(L)
        psi[r] += 1

    return psi

'''------------------------------------------------------------------------'''

def vmc(x,v):
    '''Applies a step of Variational Monte Carlo'''
    
    # Length of initial configuration
    L = np.size(x)  
    
    # Randomly select and move particle from site l to site k
    l = np.random.randint(L)
    k = np.random.randint(L)
    
        # Particle numbers on sites l and k before updates (will need later for fast computation of ratio of overlaps)
    n_l = x[l]
    n_k = x[k] 
    
    x[l] -= 1
    x[k] += 1
    
    # Reject the update if site l was vacant.
    if x[l] == -1:
        x[l] += 1
        x[k] -= 1
        return x

    #Metropolis sampling
    #REFERENCE: Becca & Sorella QMC Approaches for Correlated Systems (2017) ; Section 5.5

    #Ratio of Overlaps of new/initial with non-interacting state
    w_non_interacting = np.sqrt(n_l/(n_k+1))

    #Ratio of Jastrow factors for new/initial states
    w_jastrow = np.exp(-v*(n_k-n_l+1))

    #Weight
    w = np.abs(w_non_interacting*w_jastrow)**2
    
    # Metropolis sampling
    if w >= 1:
        # Accept update
        return x
    else:
        if np.random.random() < w:
            # Accept update
            return x
        else:
            # Reject update
            x[l] += 1
            x[k] -= 1
            return x
'''------------------------------------------------------------------------'''

def optimize_jastrow(): # This may be part of the vmc function instead
    return 0

'''------------------------------------------------------------------------'''

def bh_kinetic(x,t,v):
    '''Give a state and apply the kinetic operator of the BH-Model
    to determine its contribution to the total kinetic energy'''

    L = np.size(x) #Number of total sites in the configuration
    kineticSum = 0 #Initialize total kinetic energy contribution

    #Loop for bdag_i*b_j
    for i in range(L):
        j = (i+1)%(L) #Neighboring site with PBC taken into account
        
        #Store particle number in i,j before acting with creation/anihilation operators
        #print("HEY!!!", x)
        n_i = x[i]
        n_j = x[j]

        #Add to energy
        kineticSum += np.sqrt(n_i+1) * np.sqrt(n_j) * np.sqrt(n_j/(n_i+1)) * np.exp(-v*(n_i-n_j+1))

    #Loop for b_i*bdag_j
    for i in range(L):
        j = (i+1)%(L) #Neighboring site with PBC taken into account
            
        #Store particle number in i,j before acting with creation/anihilation operators
        n_i = x[i]
        n_j = x[j]

        #Add to energy
        kineticSum += np.sqrt(n_i)*np.sqrt(n_j+1) * np.sqrt(n_i/(n_j+1)) * np.exp(-v*(n_j-n_i+1))

    return -kineticSum

'''------------------------------------------------------------------------'''

def bh_potential(x,U,mu):

    L = np.shape(x)[0] #Number of total sites in the configuration

    #Interaction part
    potentialEnergy = 0
    for i in range(L):
        n_i = x[i] #particle number on i_th site
        
        potentialEnergy += (U/2)*(n_i*(n_i-1)) - mu*n_i

    return potentialEnergy

# To do:

    # Store the number of bins (rows: mc steps (one per L update attempts) , columns: energy?)
    # Need to set a bin size
    # FIND THE LOCAL MINIMA FOR VARIOUS U VALUES 
'''------------------------------------------------------------------------'''

def main():

    #Variational parameters
    #v = 0.7991
    
    #Bose-Hubbard parameters
    L = 4
    N = 4
    t = 1.0
    U = 0                      # v = 0.0000, V = 1.0000
    #U = 0.05000000000000       # v = 0.0145, V = 1.0091
    #U = 0.51164649614037       # v = 0.1121, V = 1.0404
    #U = 5.11646496140377       # v = 0.7991, V = 0.6819
    #U = 51.16464961403775      # v = 3.2310, V = 0.085
    #U = 199.05358527674870     #v = 4.6000,  V = 0.0081
    
    U = 2.25
    mu = 0

    #Observables
    energy = 0

    #Initialize a random BoseHubbard configuration
    x = random_boson_config(L,N)
    x = np.array([1,1,1,1])


    #print(x_old)

   # M = 10
    #Apply M steps of VMC to the configuration
   # for m in range(M):
   #     energy += (bh_kinetic(x,t,v) + bh_potential(x,U,mu))
        #print(energy)
   #     if m%2000==0:
   #         print(energy/(m+1))
   #     x = vmc(x,v)
        
    # Write ground state energies to disk as a function of variational parameter
    M = 11000
    #Apply M steps of VMC to the configuration
    
    Egs = [] #Store ground state energies
    vList = np.linspace(0,1.5,400) #List of variational params to be tested
    for v in vList:
        sweeps = 0                   #Count the number of times the loop is entered
        energy = 0
        for m in range(M):
            x = vmc(x,v)
            if m > 0.1*M: #Only measure statistics after equilibration
                sweeps+=1
                energy += (bh_kinetic(x,t,v) + bh_potential(x,U,mu))
                #print(energy/(sweeps))
                x = vmc(x,v)
        Egs.append(energy/sweeps)
    
    #Format the data file
    columns = np.array((vList,Egs))
    columns = columns.T
    with open("EGS_%d_%d_%.4f.dat"%(L,N,U),"w+") as data:
        np.savetxt(data,columns,delimiter="   ",fmt="%.4f",header="v      Egs     L=%d,N=%d,U=%.14f"%(L,N,U))
    

if __name__ == "__main__":
    main()
