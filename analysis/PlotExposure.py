#!/usr/bin/python
"""
PlotExposure.py

Plot discrimination significance as a function of

BJK 23/06/2017
"""

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
        'size'   : 16}
#Edit to 16 here!

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

#---------------------------------



#----Run Parameters---
mx = 50
R0 = sys.argv[1]

print " Plotting discrimination significance for lambda_n/lambda_p =  " + str(R0)+ " ..."

#---Functions----

#Read in a list of significances for a given point (pID) in parameter space
#for a given reconstruction (reconID)
def getSigvals_exp(reconID, exposure):
    #Filename for results file
    fname = "../results/" + reconID + "_exposure/Results_R=" + R0 + "_" + str(exposure) +'kg.txt'
    
    #Check if the file exists (and clean up the data if necessary)
    if (os.path.exists(fname)):
        data=np.loadtxt(fname)
        #Set infinite values to a large number...(say 10 sigma)
        data[data == float('+inf')] = 10.0
        data = data[~np.isnan(data)]
        if (len(data) == 0):
            data = [0]
    else:
        print " Error: File not found - " + fname
        data = np.zeros(1)

    return np.sort(data)

#Calculate significance (median, mean, upper, lower) for a given point and reconstruction
def getSignificance_exp(reconID, exposure, kind="Median"):
    sigvals = getSigvals_exp(reconID, exposure)
    Nsamps = len(sigvals)

    if (kind == "Mean"):
        return np.mean(sigvals)
    if (kind == "Median"):
        return np.median(sigvals)
    if (kind == "Upper"):
        return np.percentile(sigvals, 84.0)
    if (kind == "Lower"):
        return np.percentile(sigvals, 16.0)
    
    ind = int(np.round(Nsamps*0.1))
    return sigvals[ind]
    
    
#----Calculations---

Nvals = 32
exprange = np.round(np.logspace(2, 5, Nvals))

#Significances for ensemble A
sig_med_A = exprange*0.0
sig_upper_A = exprange*0.0
sig_lower_A = exprange*0.0

#Significances for ensemble D
sig_med_D = exprange*0.0 + 1e-45
sig_upper_D = exprange*0.0 + 1e-45
sig_lower_D = exprange*0.0 + 1e-45


reconID_A = "AP_ExptA_" + str(mx)
for i in range(Nvals):
    sig_med_A[i] = getSignificance_exp(reconID_A, int(exprange[i]), kind="Median") + 1e-45
    sig_upper_A[i] = getSignificance_exp(reconID_A, int(exprange[i]), kind="Upper") + 1e-45
    sig_lower_A[i] = getSignificance_exp(reconID_A, int(exprange[i]), kind="Lower") + 1e-45

reconID_D = "AP_ExptD_" + str(mx)
for i in range(Nvals):
    sig_med_D[i] = getSignificance_exp(reconID_D, int(exprange[i]), kind="Median") + 1e-45
    sig_upper_D[i] = getSignificance_exp(reconID_D, int(exprange[i]), kind="Upper") + 1e-45
    sig_lower_D[i] = getSignificance_exp(reconID_D, int(exprange[i]), kind="Lower") + 1e-45


#------Plotting---------

fig = pl.figure(figsize=(8,6))

ax1 = fig.add_subplot(111)

#Plot significance as a function of exposure (in ton years)
lmed_Si, = ax1.loglog(exprange/1e3, sig_med_A, linewidth=1.5, color='DarkGreen')
lband_Si = ax1.fill_between(exprange/1e3, sig_lower_A, sig_upper_A, color='ForestGreen', alpha = 0.4)

lmed_D, = ax1.loglog(exprange/1e3, sig_med_D, linewidth=1.5, color='DarkBlue')
lband_D = ax1.fill_between(exprange/1e3, sig_lower_D, sig_upper_D, color='Navy', alpha = 0.4)


leg = ax1.legend([(lmed_Si, lband_Si),(lmed_D, lband_D)], [r"Si", r"Ge + CaWO$_4$"], loc="lower right",frameon=False,fontsize=16.0)

#Sort out labels and ticks
ax1.text(0.12, 8.4, r"Fixed Xe + Ar exposure", ha="left")
ax1.text(0.12, 6.9, r"$m_\chi = " + str(mx) + "\,\,\mathrm{GeV}$; $f = -0.995$", ha="left")
ax1.text(0.12, 5.7, r"$\lambda_n/\lambda_p = "+ R0 + "$", ha="left")

#ax1.set_xlim(-1.00, -0.94)
ax1.set_ylim(0.5, 10)

ax1.set_yticks([0.5, 1, 2, 3, 4, 5, 10])
ax1.set_yticklabels([r'$0.5\sigma$', r'$1\sigma$',\
                r'$2\sigma$',r'$3\sigma$',r'$4\sigma$', r'$5\sigma$', r'$10\sigma$'])

ax1.set_xticks([0.1, 1, 10, 100])
ax1.set_xticklabels([0.1, 1, 10, 100])

ax1.set_ylabel('Discrimination significance',fontsize=18)
ax1.set_xlabel('Exposure [ton yr]',fontsize=18)

#Add some lines for 3, 5 sigma and 3 ton years
ax1.axvline(3, ymin=0, ymax = 1.00,linestyle='--', color='k')
ax1.axhline(3, linestyle=':', color='k',linewidth=1.5)
ax1.axhline(5, linestyle=':', color='k', linewidth=1.5)


pl.tight_layout()
pl.savefig("../plots/Exposure_R="+R0+".pdf", bbox_inches="tight")
#pl.show()