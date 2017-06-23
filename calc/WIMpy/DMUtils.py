# coding: utf-8

import numpy as np
from numpy import cos, sin
from scipy.integrate import trapz, cumtrapz, quad
from scipy.interpolate import interp1d
from numpy.random import rand
from scipy.special import sph_jn, erf

from numpy import pi


#---------------------------------------------------------

#Calculate eta, the velocity integral
def calcEta(v,vlag, sigmav,vesc):
    aplus = np.minimum((v+vlag), vesc)/(np.sqrt(2)*sigmav)
    aminus = np.minimum((v-vlag), vesc)/(np.sqrt(2)*sigmav)
    aesc = vesc/(np.sqrt(2)*sigmav)
    #return (np.pi/vlag)*(erf(aplus) - erf(aminus))
    
    vel_integral = 0
    
    N = 1.0/(erf(aesc) - np.sqrt(2.0/np.pi)*(vesc/sigmav)*np.exp(-0.5*(vesc/sigmav)**2))
    
    vel_integral = (0.5/vlag)*(erf(aplus) - erf(aminus))
    vel_integral -= (1.0/(np.sqrt(np.pi)*vlag))*(aplus - aminus)*np.exp(-0.5*(vesc/sigmav)**2)
    
    return 2*np.pi*N*vel_integral

# Calculate standard SI helm form factor
def calcSIFormFactor(E, m_N, old=False):
    #Helm

    #Define conversion factor from amu-->keV
    amu = 931.5*1e3

    #Convert recoil energy to momentum transfer q in keV
    q1 = np.sqrt(2*m_N*amu*E)

    #Convert q into fm^-1
    q2 = q1*(1e-12/1.97e-7)
    
    #Calculate nuclear parameters
    s = 0.9
    a = 0.52
    c = 1.23*(m_N**(1.0/3.0)) - 0.60
    R1 = np.sqrt(c*c + 7*pi*pi*a*a/3.0 - 5*s*s)
    
    if (old):
        R1 = np.sqrt((1.2**2)*m_N**(2.0/3.0) - 5)
 
    x = q2*R1
    J1 = np.sin(x)/x**2 - np.cos(x)/x
    F = 3*J1/x
    return (F**2)*(np.exp(-(q2*s)**2))

#Calculate differential rate as a function of couplings l
def dRdE(E, N_p, N_n, mx, l): 
  A = N_p + N_n   
  sig = (1.973e-14*1.973e-14)*4.0*(reduced_m(1.0, mx))**2.0/np.pi

  int_factor = 0.5*sig*calcSIFormFactor(E, A)*\
      ((l[0]*N_p + l[1]*N_n)**2.0 + (l[2]*N_p + l[3]*N_n)**2.0)
    
  return rate_prefactor(A, mx)*int_factor*calcEta(vmin(E, A, mx), vlag=232, sigmav=156, vesc=544)

#Calculate total number of events
def Nevents(E_min, E_max, N_p, N_n, mx, l):
    integ = lambda x: dRdE(x, N_p, N_n, mx, l)
    return quad(integ, E_min, E_max, epsrel=1e-4)[0]

def rate_prefactor(A, m_x):
    rho0 = 0.3
    mu = 1.78e-27*reduced_m(1.0, m_x)
    return 1.38413e-12*rho0/(4.0*np.pi*m_x*mu*mu)
    
#Reduced mass
def reduced_m(A, m_x):
    m_A = 0.9315*A
    return (m_A * m_x)/(m_A + m_x)
    
#Minimum velocity required for recoil of energy E
def vmin(E, m_N, m_x):
    res = E*0.0
    m_N2 = m_N*0.9315
    mu = (m_N2*m_x)/(m_N2+m_x)
    res =  3e5*np.sqrt((E/1e6)*(m_N2)/(2*mu*mu))
    return res
    