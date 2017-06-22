import numpy as np
from numpy import pi
import matplotlib.pyplot as pl
import matplotlib as mpl
from scipy.stats import chi2, norm
from scipy.integrate import quad
from scipy.special import iv,ive, erf
import os.path
import time
import matplotlib.colors as colors
import sys
from matplotlib.ticker import AutoMinorLocator, MultipleLocator,FormatStrFormatter
from scipy import ndimage
from scipy.interpolate import griddata

from tqdm import * #Progress Bar

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

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap




#----Run Parameters---
mx = 50
R0 = sys.argv[1]


#---Functions----

def getf(pID):
    vals = np.unravel_index(pID-1, (Np,Np), order='C')
    return frange[vals[0]]
    
def getR_PA(pID):
    vals = np.unravel_index(pID-1, (Np,Np), order='C')
    f = frange[vals[0]]
    return f/np.sqrt(1-f**2)

def getR(pID):
    vals = np.unravel_index(pID-1, (Np,Np), order='C')
    return Rrange[vals[1]]

def getSigvals_exp(expt, exp):
    reconID = "AP_Expt" + expt + "_" + str(mx)
    data=np.loadtxt("../results/" + reconID + "_exposure/Results_R=" + R0 + "_" + str(exp) +'kg.txt')
    #Set infinite values to a large number...(say 10 sigma)
    data[data == float('+inf')] = 10.0
    data = data[~np.isnan(data)]
    if (len(data) == 0):
        data = [0]
    return np.sort(data)
    #data = np.genfromtxt("../results/" + reconID + "/Results_p" + str(pID) +'.txt',\
    #        missing_values=" ", filling_values=np.nan)
    #return  np.sort(data[:,4])

def getSignificance_exp(expt, exp, kind="Mean"):
    sigvals = getSigvals_exp(expt, exp)
    Nsamps = len(sigvals)
    #if (Nsamps < 20):
    #    print " Point", pID, ": ", (Nsamps -20), " missing points"
    
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


def getAverageSamps():
    Nsamps = 0
    for pID in range(Np*Np):
     
        sigvals = getSigvals(pID+1)
        #print pID+1, len(sigvals)
        #print sigvals
        Nsamps += len(sigvals)
    return Nsamps*1.0/(Np*Np)
    
#----Calculations---

Nvals = 32

exprange = np.round(np.logspace(2, 5, Nvals))
sig_med_A = exprange*0.0
sig_upper_A = exprange*0.0
sig_lower_A = exprange*0.0

sig_med_D = exprange*0.0 + 1e-45
sig_upper_D = exprange*0.0 + 1e-45
sig_lower_D = exprange*0.0 + 1e-45

for i in range(Nvals):
    sig_med_A[i] = getSignificance_exp("A", int(exprange[i]), kind="Median") + 1e-45
    sig_upper_A[i] = getSignificance_exp("A", int(exprange[i]), kind="Upper") + 1e-45
    sig_lower_A[i] = getSignificance_exp("A", int(exprange[i]), kind="Lower") + 1e-45

for i in range(Nvals):
    sig_med_D[i] = getSignificance_exp("D", int(exprange[i]), kind="Median") + 1e-45
    sig_upper_D[i] = getSignificance_exp("D", int(exprange[i]), kind="Upper") + 1e-45
    sig_lower_D[i] = getSignificance_exp("D", int(exprange[i]), kind="Lower") + 1e-45

fig = pl.figure(figsize=(8,6))

ax1 = fig.add_subplot(111)

lmed_Si, = ax1.loglog(exprange/1e3, sig_med_A, linewidth=1.5, color='DarkGreen')
lband_Si = ax1.fill_between(exprange/1e3, sig_lower_A, sig_upper_A, color='ForestGreen', alpha = 0.4)

lmed_D, = ax1.loglog(exprange/1e3, sig_med_D, linewidth=1.5, color='DarkBlue')
lband_D = ax1.fill_between(exprange/1e3, sig_lower_D, sig_upper_D, color='Navy', alpha = 0.4)


leg = ax1.legend([(lmed_Si, lband_Si),(lmed_D, lband_D)], [r"Si", r"Ge + CaWO$_4$"], loc="lower right",frameon=False,fontsize=16.0)


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

ax1.axvline(3, ymin=0, ymax = 1.00,linestyle='--', color='k')

ax1.axhline(3, linestyle=':', color='k',linewidth=1.5)
ax1.axhline(5, linestyle=':', color='k', linewidth=1.5)

#ax1.text(0.2, 3.1, r"$3\sigma$", ha='left', va='bottom', fontsize=16)
#ax1.text(0.2, 5.1, r"$5\sigma$", ha='left', va='bottom', fontsize=16)

pl.tight_layout()

pl.savefig("../plots/Exposure_R="+R0+".pdf", bbox_inches="tight")

pl.show()





#ax1.text(-0.942, 1.06, r"Ensemble " + expt + "; $m_\chi = " + str(mx) + "\,\,\mathrm{GeV}$", ha="right")
#ax1.text(-0.9976, 0.62, r"Max. significance: $" + "{0:.2f}".format(np.max(sigvals)) + "\sigma$")

#ax1.set_ylabel(r'$\sqrt{\lambda_p^{D \,2} + \lambda_p^{\overline{D}\, 2}}/\sqrt{\lambda_n^{D \,2} + \lambda_n^{\overline{D}\, 2}}$')
#ax1.set_ylabel(r'$\lambda_n/\lambda_p$', fontsize=18)
#ax1.set_xlabel(r'$f =  (\lambda_p^D \lambda_n^{D} + \lambda_p^{\overline{D}} \lambda_n^{\overline{D}})/ \sqrt{(\lambda_p^{D \,2} + \lambda_p^{\overline{D}\, 2})(\lambda_n^{D \,2} + \lambda_n^{\overline{D}\, 2})}$')

#ax2.set_xlabel(r'$|\lambda_n^D/\lambda_n^{\overline{D}}|$')


#pl.figure()

#for i in range(32):
#    print i+1, "; l_n/l_p = ", getR(i+1)
#    print getSigvals(i+1)

#print pA
#sigs = getSigvals(pA)
#print sigs
#pl.hist(sigs)


#pl.savefig("../plots/Exposure-"+reconID + ".pdf", bbox_inches="tight")
#pl.show()