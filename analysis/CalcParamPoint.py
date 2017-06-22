import numpy as np
from numpy import pi

#Code for calculating couplings (given values of Rnp and f, or an index)

#Define the grid of coupling values to evaluate over
Np = 32
Rrange = np.linspace(0.6, 1.00, Np)
frange = np.append(np.linspace(-1.0, -0.97, 3*Np/4 + 1)[:-1], np.linspace(-0.97, -0.94,Np/4))

#Format of coupling vector is always:
#l_p_D, l_n_D, l_p_Dbar, l_n_Dbar                                                                                                                             

#Calculate cp, cn and ratios from list of couplings 'lam'
                                                                                                                                         
def Calc_cp(lam):
    return np.sqrt((lam[0]**2 + lam[2]**2)/2.0)

def Calc_cn(lam):
    return np.sqrt((lam[1]**2 + lam[3]**2)/2.0)

def Calc_Rnp(lam):
    return Calc_cn(lam)/Calc_cp(lam)
                                                                                                                                       
def Calc_f(lam):
    return (lam[0]*lam[1] + lam[2]*lam[3])*0.5/(Calc_cp(lam)*Calc_cn(lam))


#Generate list of couplings from an index or values of (Rnp, f)                                                                                      
def GeneratePoint_ind(index):
    vals = np.unravel_index(index-1, (Np,Np), order='C')
    fval = frange[vals[0]]
    Rval = Rrange[vals[1]]
    
    return GeneratePoint(Rval, fval)
                                                  
def GeneratePoint(Rnp, f):
    #Assume that l_p_D = 1 and l_p_Dbar = 0                                                                                                                   
    l_p_D = 1.0
    l_p_Dbar = 0.0

    l_n_D = f*Rnp
    l_n_Dbar = np.sqrt((1.0-f**2))*Rnp

    return np.array([l_p_D, l_n_D, l_p_Dbar, l_n_Dbar])

#Get values of f and R from a linear index

def getf(index):
    vals = np.unravel_index(index-1, (Np,Np), order='C')
    return frange[vals[0]]
    
#Particle-Antiparticle ratio
def getR_PA(index):
    vals = np.unravel_index(index-1, (Np,Np), order='C')
    f = frange[vals[0]]
    return f/np.sqrt(1-f**2)

#Proton-Neutron ratio
def getR(index):
    vals = np.unravel_index(index-1, (Np,Np), order='C')
    return Rrange[vals[1]]
