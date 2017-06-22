#!/usr/bin/python

# PlotContours_row.py
#
# Plot contours of discrimination significance over a range of DM masses
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



#----Run Parameters---

if (len(sys.argv) != 2):
    print " PlotContours_row.py requires 1 argument - e.g. PlotContours_row.py EXPT"
    print "     where EXPT = A, B, C, D"
    print " Exiting..."
    sys.exit()
    
expt = str(sys.argv[1])

print " Plotting row of contour-plots for ensemble " + expt + "..."

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

#Read in a list of significances for a given point (pID) in parameter space
#for a given reconstruction (reconID)
def getSigvals(reconID, pID):
    #Filename for results file
    fname = "../results/" + reconID + "/Results_p" + str(pID) +'.txt'
    
    #Check if the file exists (and clean up the data if necessary)
    if (os.path.exists(fname)):
        data=np.loadtxt(fname)
        data = data[data != float('+inf')]
        data = data[~np.isnan(data)]
        if (len(data) == 0):
            data = [0]
    else:
        print " Error: File not found - " + fname
        data = np.zeros(1)

    return np.sort(data)

#Calculate significance (median, mean, upper, lower) for a given point and reconstruction
def getSignificance(reconID, pID, kind="Median"):
    sigvals = getSigvals(reconID, pID)
    Nsamps = len(sigvals)

    if (kind == "Mean"):
        return np.mean(sigvals)
    if (kind == "Median"):
        return np.median(sigvals)
    if (kind == "Upper"):
        ind = int(np.round(Nsamps*(1.0-0.32)))
        return sigvals[ind]
    if (kind == "Lower"):
        ind = int(np.round(Nsamps*0.32))
        return sigvals[ind]
    
    ind = int(np.round(Nsamps*0.1))
    return sigvals[ind]

#Calculating top axis from bottom axis
def inv_tick_function(X):
    return X*1.0/np.sqrt(1+X**2)


#----Calculations---


fig, axarr = pl.subplots(1,3, figsize=(18,6))

levels = np.array([1,2,3,4,5]) #Contours of sigma
colmap = pl.get_cmap("Greens") #Color map
masses = [50, 300, 1000] #Plot for three masses

for j in range(3):

    mx = masses[j]
    reconID = "AP_Expt" + expt + "_" + str(mx)
    print "    m_x = " + str(mx) + " GeV..."

    #Number of grid points in the parameter space in each direction
    #CPP (i.e. CalcParamPoint.py) transforms indices into values of couplings
    Np = CPP.Np
    xvals = CPP.frange #x-axis
    yvals = CPP.Rrange #y-axis
    
    sigvals = np.zeros(Np*Np)
    #Get significance for each point
    for i in range(Np*Np):
        sigvals[i] = getSignificance(reconID, i+1, kind="Median")
        
    zvals = np.reshape(sigvals, (Np, Np)).T

    #Resample and filter the median results along the y-axis
    Nvals = 50
    xgrid = np.linspace(-1.0, -0.94, Nvals)
    ygrid = np.linspace(0.6, 1.00, Nvals)
    z_interp = interp2d(xvals, yvals, zvals, kind='linear')

    xvals = xgrid*1.0
    yvals = ygrid*1.0
    zvals = z_interp(xvals, yvals)

    for i in range(Nvals):
        zvals[:,i] = ndimage.filters.median_filter(zvals[:,i], 5)

    # Do some plotting
    ax1 = axarr.flatten()[j]
    
    #Plot filled contours
    cf = ax1.contourf(xvals, yvals, zvals, \
        levels, cmap=colmap, extend='max')
    #Plot contour lines
    cons0 = ax1.contour(xvals, yvals, zvals, \
        levels, colors='forestgreen')

    #Find and plot maximum point
    maxID = np.argmax(sigvals)
    ax1.plot(CPP.getf(maxID+1), CPP.getR(maxID+1),ms=12,marker='*', mew=0.5, color='k')
    print "     Maximum significance: ", np.max(sigvals), "[INDEX " + str(maxID+1)+"]"

    #Add red squares in some cases
    if ((expt == "D") and (mx == 50)):
        ax1.plot(-0.995, 0.75, ms=8,marker='s',color='r', mew=0.5)
        ax1.plot(-0.995, 0.8, ms=8,marker='s',color='r', mew=0.5)

    ax1.set_xlim(-1.00, -0.94)
    ax1.set_ylim(0.5, 1.1)
    ax1.yaxis.set_major_locator(MultipleLocator(0.1))

    #Add horizontal dashed lines for different elements
    ax1.axhline(1.0/R_Xe, linestyle="--", color='k')
    ax1.text(-0.943, 1.0/R_Xe+0.008, r"Xe",ha="right")
    ax1.axhline(1.0/R_Ar, linestyle="--", color='k')
    ax1.text(-0.943, 1.0/R_Ar+0.008, r"Ar",ha="right")

    if (expt == "A"):
        ax1.axhline(1.0/R_Si, linestyle="--", color='k')
        ax1.text(-0.943, 1.0/R_Si-0.03, r"Si",ha="right")

    if (expt == "B" or expt == "D"):
        ax1.axhline(1.0/R_Ge, linestyle="--", color='k')
        ax1.text(-0.943, 1.0/R_Ge+0.008, r"Ge",ha="right")
    
    if (expt == "C" or expt == "D"):
        ax1.axhline(1.0/R_Ca, linestyle="--", color='k')
        ax1.text(-0.943, 1.0/R_Ca-0.03, r"Ca,O", ha="right")
        ax1.axhline(1.0/R_W, linestyle="--", color='k')
        ax1.text(-0.943, 1.0/R_W-0.03, r"W", ha="right")

    #Sort out the top axis
    ax2 = ax1.twiny()
    topticks = np.array([-30,-10, -5,-4,-3])
    ax2.set_xticks(inv_tick_function(topticks))
    ax2.set_xticklabels(np.abs(topticks))
    ax2.set_xlim(ax1.get_xlim())
    
    
    #Add some labels
    ax1.text(-0.942, 1.055, r"Ensemble " + expt + "; $m_\chi = " + str(mx) + "\,\,\mathrm{GeV}$", ha="right", fontsize=18.0)
    ax1.text(-0.999, 0.52, r"Max. significance ($\bigstar$): $" + "{0:.1f}".format(np.max(sigvals)) + "\sigma$", fontsize=16.0)

    ax1.set_ylabel(r'$\lambda_n/\lambda_p$', fontsize=20)
    ax2.set_xlabel(r'$|\lambda_n^D/\lambda_n^{\overline{D}}|$', fontsize=20.0)
    
    #Deal with the overlapping subplot axis labels
    if (j in [1,2]):
        ax1.set_yticklabels([])
        ax1.set_ylabel("")
    if (j < 2):
        ax1.set_xticklabels(["-1", "-0.99", "-0.98", "-0.97", "-0.96","-0.95"])
    if (j == 2):
        ax1.set_xticklabels(["-1", "-0.99", "-0.98", "-0.97", "-0.96","-0.95", "-0.94"])


fig.suptitle( r'$f =  (\lambda_p^D \lambda_n^{D} + \lambda_p^{\overline{D}} \lambda_n^{\overline{D}})/ \sqrt{(\lambda_p^{D \,2} + \lambda_p^{\overline{D}\, 2})(\lambda_n^{D \,2} + \lambda_n^{\overline{D}\, 2})}$', \
        ha='center',x=0.5, y=0.05, fontsize=20.0)

#Add colorbar
cbar_ax = fig.add_axes([0.96, 0.15, 0.015, 0.7])
cb0 = fig.colorbar(cf, cax=cbar_ax, ticks=levels, extend='max')
cb0.set_ticklabels([r'$1\sigma$',\
        r'$2\sigma$',r'$3\sigma$',r'$4\sigma$',r'$5\sigma$'])
cb0.ax.tick_params(labelsize=18.0)

#Adjust and save to file
fig.subplots_adjust(hspace=0.1, wspace=0.05, right=0.95)
pl.savefig("../plots/Contours_row-" + expt+ ".pdf", bbox_inches="tight")
pl.show()