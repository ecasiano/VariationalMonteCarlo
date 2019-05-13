#Build Hamiltonian Matrix representing the Bose-Hubbard Model

from sympy.utilities.iterables import multiset_permutations,ordered_partitions
import numpy as np
import argparse

'''---Command Line Arguments---'''
#Create an ArgumentParser object
parser = argparse.ArgumentParser()

#Add Arguments
parser.add_argument("--bosons",action="store_true")
parser.add_argument("--fermions",action="store_true")
args = parser.parse_args()

'''---Define Functions'''

def bosonicConfigurations(L,N):
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

def fermionicConfigurations(L,N):
    '''Input: 1D Lattice Size and Number of Fermions
    Output: All possible configurations of Fermions'''
    
    configurations = [0]
    configurations = L*[0]
    
    for i in range(N):
        configurations[i] = 1

    configurations = list(multiset_permutations(configurations))
    
    return np.array(configurations)

def boseHubbardKinetic(bra,ket,t=1):
    '''Give a state and apply the kinetic operator of the BH-Model
    to determine its contribution to the total kinetic energy'''

    m = np.copy(ket)
    L = np.shape(ket)[0] #Number of total sites in the configuration
    kineticSum = 0 #Initialize total kinetic energy contribution
    
    #Loop for bdag_i*b_{i+1}
    for i in range(L):
        #Create boson on site i.
        m[i] += 1

        #Annihilate boson on site i+1 (do nothing if no bosons)
        if m[(i+1)%(L)] != 0: m[(i+1)%(L)] -= 1 
        else: m[(i+1)%(L)] = 0
        
        if np.array_equal(bra,m):
            kineticSum += np.sqrt(ket[i]+1)*np.sqrt(ket[(i+1)%(L)])
            #print("EQUAL")
            #print(ket,m)
        #else:
            #print("DIFFERENT")
            #print(ket,m)
        
        #Make m a copy of the original state n again.
        m = np.copy(ket)
    
    #Loop for b_i*bdag_{i+1}
    for i in range(L):
        #Create boson on site i+1.
        m[(i+1)%(L)] += 1

        #Annihilate boson on site i (do nothing if no bosons)
        if m[i] != 0: m[i] -= 1 
        else: m[i] = 0

        if np.array_equal(bra,m):
            kineticSum += np.sqrt(ket[i])*np.sqrt(ket[(i+1)%(L)]+1)
            #print("EQUAL")
            #print(ket,m)
        #else:
            #print("DIFFERENT")
            #print(ket,m)

        #Make m a copy of the original state n again.
        m = np.copy(ket)
        
    return -kineticSum

def tVKinetic(bra,ket,t=1):
    '''Give a state and apply the kinetic operator of the BH-Model
    to determine its contribution to the total kinetic energy'''

    m = np.copy(ket)
    L = np.shape(ket)[0] #Number of total sites in the configuration
    kineticSum = 0 #Initialize total kinetic energy contribution
    
    #Loop for bdag_i*b_{i+1}
    for i in range(L):
        #Create boson on site i.
        if m[i] == 0:
            m[i] += 1
        else: m[i] = 0

        #Annihilate boson on site i+1 (do nothing if no bosons)
        if m[(i+1)%(L)] != 0: m[(i+1)%(L)] -= 1 
        else: m[(i+1)%(L)] = 0
        
        if np.array_equal(bra,m):
            kineticSum += np.sqrt(ket[i]+1)*np.sqrt(ket[(i+1)%(L)])
            #print("EQUAL")
            #print(ket,m)
        #else:
            #print("DIFFERENT")
            #print(ket,m)
        
        #Make m a copy of the original state n again.
        m = np.copy(ket)
    
    #Loop for b_i*bdag_{i+1}
    for i in range(L):
        #Create boson on site i+1.
        if m[(i+1)%(L)] == 0:
            m[(i+1)%(L)] += 1
        else: m[(i+1)%(L)] = 0

        #Annihilate boson on site i (do nothing if no bosons)
        if m[i] != 0: m[i] -= 1 
        else: m[i] = 0

        if np.array_equal(bra,m):
            kineticSum += 1
            #print("EQUAL")
            #print(ket,m)
        #else:
            #print("DIFFERENT")
            #print(ket,m)

        #Make m a copy of the original state n again.
        m = np.copy(ket)
        
    return -kineticSum

def boseHubbardHamiltonian(configurations):
    '''Input: Set of all possible configurations of bosons on a 1D Lattice'''
    
    #Store HilbertSpace Size
    hilbertSize = np.shape(configurations)[0]
    #print(hilbertSize)
    
    #Store Lattice Size
    L = np.shape(configurations)[1]
    
    #Initialize Hamiltonian Matrix
    H = np.zeros((hilbertSize,hilbertSize))
    
    #Print out configurations to determine ordering of the basis
    #print("Configurations: ", configurations)
    #Fill in upper diagonal of the Hamiltonian
    for i in range(hilbertSize):
        bra = configurations[i]
        for j in range(i,hilbertSize):
            ket = configurations[j]
            H[i,j] = boseHubbardKinetic(bra,ket)
            H[j,i] = H[i,j] #Use Hermiticity to fill up lower diagonal
            
    return H

def tVHamiltonian(configurations):
    '''Input: Set of all possible configurations of bosons on a 1D Lattice'''
    
    #Store HilbertSpace Size
    hilbertSize = np.shape(configurations)[0]
    #print(hilbertSize)
    
    #Store Lattice Size
    L = np.shape(configurations)[1]
    
    #Initialize Hamiltonian Matrix
    H = np.zeros((hilbertSize,hilbertSize))
    
    #Print out configurations to determine ordering of the basis
    #print("Configurations: ", configurations)
    #Fill in upper diagonal of the Hamiltonian
    for i in range(hilbertSize):
        bra = configurations[i]
        for j in range(i,hilbertSize):
            ket = configurations[j]
            H[i,j] = tVKinetic(bra,ket)
            H[j,i] = H[i,j] #Use Hermiticity to fill up lower diagonal
            
    return H

def Pn_BH(psi,configurations,lA):
    '''Input: Quantum State in second quantization bipartitioned into 
    subregions of size lA and L-lA'''
    '''Output: Probabilities of measuring states with particle number n in subregion A'''  
    L = np.shape(configurations[-1])[0]
    N = configurations[-1,0]
    hilbertSize = np.shape(configurations)[0]

    #Array to store probabilities
    Pn = np.zeros(N+1)
    for n in range(N+1):
        for i in range(hilbertSize):
            if np.sum(configurations[i][0:lA]) == n:
                Pn[n] += psi[i]**2

    return Pn

def Pn_tV(psi,configurations,lA):
    '''Input: Quantum State in second quantization bipartitioned into 
    subregions of size lA and L-lA'''
    '''Output: Probabilities of measuring states with particle number n in subregion A'''  
    L = np.shape(configurations[-1])[0]
    N = np.sum(configurations[-1])
    hilbertSize = np.shape(configurations)[0]

    #Array to store probabilities
    Pn = np.zeros(N+1)
    for n in range(N+1):
        for i in range(hilbertSize):
            if np.sum(configurations[i][0:lA]) == n:
                Pn[n] += psi[i]**2

    return Pn
       
'''---Main---'''
def main():
    
    #BUG REPORT: ONLY WORKS FOR N <= L at the moment
    #Parameters
    L = 6
    N = 3
    lA = 3
    t = 1
    V = 0
    
    '''---Bosons---'''
    if args.bosons:
        #Store all possible configurations of N bosons on L lattice sites
        configurations = bosonicConfigurations(L,N)
        print("Configurations: ",configurations)

        #Hamiltonian
        H = (boseHubbardHamiltonian(configurations))

        print("")
        print("L=%d, N=%d, t=1, V = 0"%(L,N))
        print("Hamiltonian:")
        print("")
        print(H)

        #Find ground state energy and state of the Hamiltonian
        eigs,evecs = np.linalg.eigh(H)
        #print(eigs)
        egs = eigs[0]
        psi = evecs[:,0]
        print("")
        print("Ground State Energy: ",egs)
        print("Ground State:   ",psi)
        print("norm(psi): ",np.linalg.norm(psi))
        print("")

        #Determine the probabilities for each particle number sector in n
        n = np.arange(N+1)
        probs = Pn_BH(psi,configurations,lA)

        #Print probabilities of getting particle number n in partition A
        print("n:  ", n)
        print("Pn: ", probs)

        #Save probabilities to file
        data = np.c_[n,probs]
        fileName = "pnL%dN%dlA%dBH.dat"%(L,N,lA)
        header = "L = %d, N = %2d, lA = %d, t = %f, V = %f (Bosons)\nn      Pn"%(L,N,lA,t,V)
        file = np.savetxt(fileName,data,fmt="%f",header=header)

        print(np.sum(probs))
        
    '''---Fermions---'''
    if args.fermions:
        #Store all possible configurations of N bosons on L lattice sites
        configurations = fermionicConfigurations(L,N)
        print("Configurations: ",configurations)

        #Hamiltonian
        H = (tVHamiltonian(configurations))

        print("")
        print("L=%d, N=%d, t=1, V = 0"%(L,N))
        print("Hamiltonian:")
        print("")
        print(H)

        #Find ground state energy and state of the Hamiltonian
        eigs,evecs = np.linalg.eigh(H)
        #print(eigs)
        egs = eigs[0]
        psi = evecs[:,0]
        print("")
        print("Ground State Energy: ",egs)
        print("Ground State:   ",psi)
        print("norm(psi): ",np.linalg.norm(psi))
        print("")

        #Determine the probabilities for each particle number sector in n
        n = np.arange(N+1)
        probs = Pn_tV(psi,configurations,lA)

        #Print probabilities of getting particle number n in partition A
        print("n:  ", n)
        print("Pn: ", probs)

        #Save probabilities to file
        data = np.c_[n,probs]
        fileName = "pnL%dN%dlA%dtV.dat"%(L,N,lA)
        header = "L = %d, N = %2d, lA = %d, t = %f, V = %f (Fermions)\nn      Pn"%(L,N,lA,t,V)
        file = np.savetxt(fileName,data,fmt="%f",header=header)

        print(np.sum(probs))

if __name__ == "__main__":
    main()   