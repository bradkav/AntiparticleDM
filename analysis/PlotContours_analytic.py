#!/usr/bin/python

# PlotContours_analytic.py
#
# Plot contours of discrepancy between Dirac and Majorana cross sections
# calculated using AnalyticEstimate.nb
#
# BJK 22/06/2017

import numpy as np
from numpy import pi
import matplotlib.pyplot as pl
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib.ticker import MultipleLocator
from scipy import ndimage

#------ Matplotlib parameters ------

import matplotlib.pyplot as pl
import matplotlib as mpl
import matplotlib.colors as colors

font = {'family' : 'sans-serif',
        'size'   : 17}

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

#-----------------------------------



#---Constants-----

# Proton-to-neutron ratios for different nuclei
R_Xe = (131.0-54)/54
R_Ar = (40-18.0)/18.0
R_Ge = (73.0-32.0)/32.0
R_Si = 1.0
R_Ca = 1.0
R_O = 1.0
R_W = (184-74.0)/74

#---Functions----
    
#Calculating top axis from bottom axis
def inv_tick_function(X):
    return X*1.0/np.sqrt(1+X**2)
    
#----Calculations---


#Import results calculated using AnalyticEstimate.nb notebook
xvals, yvals, zvals = np.loadtxt('../results/AnalyticEstimate.dat', unpack=True)

xvals = np.reshape(xvals, (61,61))
yvals = np.reshape(yvals, (61,61))
zvals = np.reshape(zvals, (61,61))

#The calculations were done for Ensemble A:
expt = "A"

colmap = pl.get_cmap("Purples")
levels = np.linspace(0.1, 1.0,6)

fig = pl.figure(figsize=(8,6))
ax1 = fig.add_subplot(111)

#Plot contours of Delta
cf = ax1.contourf(xvals, yvals, zvals, \
    levels, cmap=colmap, extend='max')
cons0 = ax1.contour(xvals, yvals, zvals, \
    levels, colors='DarkMagenta')
cb0 = pl.colorbar(cf, ticks=levels, extend='max')

ax2 = ax1.twiny()

ax1.axhline(1.0/R_Xe, linestyle="--", color='k')
ax1.text(-0.943, 1.0/R_Xe+0.008, r"Xe",ha="right")
ax1.axhline(1.0/R_Ar, linestyle="--", color='k')
ax1.text(-0.943, 1.0/R_Ar+0.008, r"Ar",ha="right")
ax1.axhline(1.0/R_Si, linestyle="--", color='k')
ax1.text(-0.943, 1.0/R_Si-0.03, r"Si",ha="right")

ax1.set_xlim(-1.00, -0.94)
ax1.set_ylim(0.5, 1.1)

ax1.yaxis.set_major_locator(MultipleLocator(0.1))

#Sort out the top axis
topticks = np.array([-30,-10, -5,-4,-3])
ax2.set_xticks(inv_tick_function(topticks))
ax2.set_xticklabels(np.abs(topticks))
ax2.set_xlim(ax1.get_xlim())


#Add some labels
ax1.text(-0.942, 1.06, r"Ensemble " + expt, ha="right")
ax1.text(-0.9976, 0.53, r"Cross section discrepancy: $(\sigma^D - \sigma^M)^2/(\sigma^D)^2$", ha="left", fontsize=15.0)

ax1.set_ylabel(r'$\lambda_n/\lambda_p$', fontsize=22.0)
ax1.set_xlabel(r'$f =  (\lambda_p^D \lambda_n^{D} + \lambda_p^{\overline{D}} \lambda_n^{\overline{D}})/ \sqrt{(\lambda_p^{D \,2} + \lambda_p^{\overline{D}\, 2})(\lambda_n^{D \,2} + \lambda_n^{\overline{D}\, 2})}$', fontsize=20)

ax2.set_xlabel(r'$|\lambda_n^D/\lambda_n^{\overline{D}}|$', fontsize=20.0)

pl.savefig("../plots/Significance-Estimate.pdf", bbox_inches="tight")
pl.show()