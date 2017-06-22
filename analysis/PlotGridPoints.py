#!/usr/bin/python

# PlotGridPoints.py
#
# Plot DM-coupling grid points as a function of point ID
#
# BJK 22/06/2017

import numpy as np
from numpy import pi
from scipy.integrate import quad
from scipy.interpolate import interp1d, interp2d
from scipy import ndimage
from matplotlib.ticker import MultipleLocator
import os.path
import sys

import CalcParamPoint as CPP

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


#---Functions-----

#Calculating top axis from bottom axis
def inv_tick_function(X):
    return X*1.0/np.sqrt(1+X**2)


#----Calculations---


fig, ax1 = pl.subplots(1,1, figsize=(7,6))

#Number of grid points in the parameter space in each direction
#CPP (i.e. CalcParamPoint.py) transforms indices into values of couplings
Np = CPP.Np
xvals = CPP.frange #x-axis
yvals = CPP.Rrange #y-axis
indices = np.arange(1, Np*Np+1)
#print indices.shape

xgrid, ygrid = np.meshgrid(xvals, yvals)
zgrid = np.reshape(indices, (Np, Np)).T

#Plot grid points
cf = ax1.scatter(xgrid, ygrid, c=zgrid)

#Sort out the top axis
ax2 = ax1.twiny()
topticks = np.array([-30,-10, -5,-4,-3])
ax2.set_xticks(inv_tick_function(topticks))
ax2.set_xticklabels(np.abs(topticks))
ax2.set_xlim(ax1.get_xlim())

ax1.set_ylabel(r'$\lambda_n/\lambda_p$', fontsize=20)
ax2.set_xlabel(r'$|\lambda_n^D/\lambda_n^{\overline{D}}|$', fontsize=20.0)

ax1.set_xlim(-1.01, -0.93)
ax1.set_xticklabels([" ","-1", "-0.99", "-0.98", "-0.97", "-0.96","-0.95", "-0.94", " "])

ax1.set_xlabel( r'$f =  (\lambda_p^D \lambda_n^{D} + \lambda_p^{\overline{D}} \lambda_n^{\overline{D}})/ \sqrt{(\lambda_p^{D \,2} + \lambda_p^{\overline{D}\, 2})(\lambda_n^{D \,2} + \lambda_n^{\overline{D}\, 2})}$', \
        ha='center',x=0.5, y=0.05, fontsize=20.0)

ax1.text(-1.002, 0.57, "1")
ax1.text(-0.94, 1.01, "1024")

#Add colorbar
tickvals = np.append([1],np.arange(128,1025, 128))
cbar_ax = fig.add_axes([0.95, 0.15, 0.03, 0.7])
cb0 = fig.colorbar(cf, cax=cbar_ax,ticks=tickvals)
cb0.set_label('Point ID')

pl.savefig("../plots/GridPoints.pdf", bbox_inches="tight")
pl.show()