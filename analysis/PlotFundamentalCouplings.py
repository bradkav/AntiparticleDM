#!/usr/bin/python
"""
PlotFundamentalCouplings.py

Plot an illustration of the fundamental (Lagrangian-level) 
couplings which correspond to the regions where good 
discrimination is possible.

BJK 15/09/2017
"""

import numpy as np

#----- Matplotlib parameters ------

import matplotlib.pyplot as pl
import matplotlib as mpl
import matplotlib.colors as colors


font = {'family' : 'sans-serif',
        'size'   : 18}
mpl.rc('font', **font)

mpl.rcParams['xtick.major.size'] = 8
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['xtick.minor.size'] = 3
mpl.rcParams['xtick.minor.width'] = 1
mpl.rcParams['ytick.major.size'] = 8
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['ytick.minor.size'] = 3
mpl.rcParams['ytick.minor.width'] = 1
mpl.rcParams['axes.linewidth'] = 1.5

mpl.rc('font', **font)
mpl.rc('text.latex', preamble=[r'\usepackage{color}', r'\usepackage{amssymb}'])

#-------------------------------

#---Functions---

#Calculate fundamental Lagrangian couplings (see Eq. 1)
#for a given value of f and rnp = lambda_n/lambda_p
#Note: lne denotes lambda_{n,e} etc. (coupling of
#even operator to neutrons).
#Returns couplings in format (lpe, lpo, lne, lno)

#Positive Solution
def CalculateFundamental_plus(rnp, f, lne, lno = 1.0):
    
    theta = np.arccos(f)
    
    lpe = (np.cos(theta)*lne - np.sin(theta)*lno)*1.0/rnp
    lpo = (np.sin(theta)*lne + np.cos(theta)*lno)*1.0/rnp

    return (lpe, lpo, lne, lno)

#Negative Solution
def CalculateFundamental_minus(rnp, f, lne, lno = 1.0):
    
    theta = -np.arccos(f)
    
    lpe = (np.cos(theta)*lne - np.sin(theta)*lno)*1.0/rnp
    lpo = (np.sin(theta)*lne + np.cos(theta)*lno)*1.0/rnp

    return (lpe, lpo, lne, lno)
    

#----Calculations----

# Generate 100 samples each with a different value of lne
# (lno is fixed to 1)
# These 100 samples will correspond to 200 fundamental couplings
# because there is a positive and negative solution for each.
Np = 100
points = np.zeros((2*Np, 4))
lne_vals = np.linspace(-1.0, 1.0, Np)

# Specify central values (and spread) for r_np and f
rnp0 = 0.75
f0 = -0.99

delta_rnp = 0.05
delta_f = 0.005

#Calculate fundamental couplings (for values of r_np and f near the central values)
for i in range(Np):
    points[2*i,:] = CalculateFundamental_plus(rnp0+delta_rnp*(2.0*np.random.rand()-1.0), f0+delta_f*(2.0*np.random.rand()-1.0), lne_vals[i])
    points[2*i+1,:] = CalculateFundamental_minus(rnp0+delta_rnp*(2.0*np.random.rand()-1.0), f0+delta_f*(2.0*np.random.rand()-1.0), lne_vals[i])
    
#----Plotting----

fig, ax = pl.subplots(1,1, figsize=(8,6))

scat = ax.scatter(points[:,0], points[:,2], c=points[:,1],vmin=-1.5, vmax=-1.0)


# Add diagonal line for central value of r_np
ax.plot([-2.0, 2.0], [2.0*0.75, -2.0*0.75], 'k--', linewidth=2.0)

# Some grid lines
ax.axvline(0, linestyle='-', color='grey', alpha=0.5)
ax.axhline(0, linestyle='-', color='grey', alpha=0.5)
ax.set_xlim(-1.75,1.75)
ax.set_ylim(-1.75,1.75)

# Add colorbar
cbar = fig.colorbar(scat,ax=ax)
cbar.ax.set_ylabel("$\lambda_{p,o}$", fontsize=20)

# Add labels
ax.set_xlabel("$\lambda_{p,e}$", fontsize=20)
ax.set_ylabel("$\lambda_{n,e}$", fontsize=20)
ax.text(-1.55,-1.0, r"$f \sim - 0.99$", ha="left")
ax.text(-1.55,-1.25, r"$\lambda_n/\lambda_p \sim 0.75$", ha="left")
ax.text(-1.55,-1.5, r"$\lambda_{n,o} = 1$", ha="left")

# Output to file
pl.tight_layout()
pl.savefig("../plots/FundamentalCouplings.pdf", bbox_inches="tight")
#pl.show()



