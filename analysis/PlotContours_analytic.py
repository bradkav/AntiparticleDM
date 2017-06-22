import numpy as np
from numpy import pi
import matplotlib.pyplot as pl
import matplotlib as mpl
from scipy.stats import chi2, norm
from scipy.integrate import quad
from scipy.special import iv,ive, erf
from scipy.interpolate import interp1d, interp2d
import os.path
import time
import matplotlib.colors as colors
import sys
from matplotlib.ticker import AutoMinorLocator, MultipleLocator,FormatStrFormatter
from scipy import ndimage
from scipy.interpolate import griddata

from tqdm import * #Progress Bar

font = {'family' : 'sans-serif',
        'size'   : 17}
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

#mpl.rc('text',usetex=True)
mpl.rc('text.latex', preamble=[r'\usepackage{color}', r'\usepackage{amssymb}'])

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap




#----Run Parameters---
 


#---Constants-----
Np = 32
#Rrange = np.linspace(0.6, 1.5, Np)
#frange = np.linspace(-1.0, -0.95, Np)                                        
#vals = np.unravel_index(pID-1, (Np,Np), order='C')                           
#frange = np.append(np.linspace(-1.0, -0.975, 3*Np/4), np.linspace(-0.975, -0\
#.95, Np/4))                                                                   
#frange = np.linspace(-1.0, -0.8, Np)

Rrange = np.linspace(0.6, 1.00, Np)
#frange = np.linspace(-1.0, -0.95, Np)                                        
#vals = np.unravel_index(pID-1, (Np,Np), order='C')                           
frange = np.append(np.linspace(-1.0, -0.97, 3*Np/4 + 1)[:-1], np.linspace(-0.97, -0.94,Np/4))

#print frange

R_Xe = (131.4-54)/54
R_Ar = (40-18.0)/18.0

R_Ge = (73.0-32.0)/32.0
R_Si = 1.0

R_Ca = 1.0
R_O = 1.0
R_W = (184-74.0)/74

R_CaWO4 = 0.139*R_Ca + 0.222*R_O + 0.639*R_W

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

def getSigvals(pID):
    
    fname = "../results/" + reconID + "/Results_p" + str(pID) +'.txt'
    
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
    #data = np.genfromtxt("../results/" + reconID + "/Results_p" + str(pID) +'.txt',\
    #        missing_values=" ", filling_values=np.nan)
    #return  np.sort(data[:,4])

def getSignificance(pID, kind="Mean"):
    sigvals = getSigvals(pID)
    Nsamps = len(sigvals)
    #if (Nsamps < 20):
    #    print " Point", pID, ": ", (Nsamps -20), " missing points"
    
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
    
def getAverageSamps():
    Nsamps = 0
    for pID in range(Np*Np):
     
        sigvals = getSigvals(pID+1)
        #print pID+1, len(sigvals)
        #print sigvals
        Nsamps += len(sigvals)
    return Nsamps*1.0/(Np*Np)
    
#----Calculations---

pA =40+32

xvals, yvals, zvals = np.loadtxt('AnalEstimate.dat', unpack=True)

xvals = np.reshape(xvals, (61,61))
yvals = np.reshape(yvals, (61,61))
zvals = np.reshape(zvals, (61,61))


#colmap = truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100)
cm = pl.get_cmap("Purples")
colmap = truncate_colormap(cm, 0.0, 1.00)

#pl.figure()
#pl.scatter(xvals, yvals, c=sigvals)
#pl.colorbar()
#pl.xlim(-1.00, -0.95)
#pl.ylim(0.75, 1.75)

#levels = np.array([1,2,3,4,5])
#levels = np.logspace(-1, np.log10(2.0),5)
levels = np.linspace(0.1, 1.0,6)

fig = pl.figure(figsize=(8,6))

ax1 = fig.add_subplot(111)


cf = ax1.contourf(xvals, yvals, zvals, \
    levels, cmap=colmap, extend='max')
cons0 = ax1.contour(xvals, yvals, zvals, \
    levels, colors='DarkMagenta')
cb0 = pl.colorbar(cf, ticks=levels, extend='max')
#cb0.set_ticklabels([r'$1\sigma$',\
#        r'$2\sigma$',r'$3\sigma$',r'$4\sigma$',r'$5\sigma$'])


#pl.clabel(cons0, cons0.levels, inline=True, \
#            fmt = r'%d $\sigma$', fontsize=12, fontcolor='k')

#ax1.plot(getf(maxID+1), getR(maxID+1),ms=12,marker='*', mew=0.5, color='k')
#ax1.plot(getf(testind+1), getR(testind+1),ms=12,marker='+', mew=0.2)

expt = "A"

if ((expt == "D") and (mx == 50)):
    ax1.plot(-0.995, 0.75, ms=8,marker='s',color='r', mew=0.5)
    ax1.plot(-0.995, 0.8, ms=8,marker='s',color='r', mew=0.5)

ax2 = ax1.twiny()

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
    #ax1.axhline(1.0/R_CaWO4, linestyle="--", color='k')
    #ax1.text(-0.942, 1.0/R_CaWO4-0.03, r"CaWO$_4$", ha="right")
    ax1.axhline(1.0/R_Ca, linestyle="--", color='k')
    ax1.text(-0.943, 1.0/R_Ca-0.03, r"Ca,O", ha="right")
    ax1.axhline(1.0/R_W, linestyle="--", color='k')
    ax1.text(-0.943, 1.0/R_W-0.03, r"W", ha="right")

#ax1.set_xlim(-20.00, -1.0)
ax1.set_xlim(-1.00, -0.94)
ax1.set_ylim(0.5, 1.1)

def inv_tick_function(X):
    #V = X*1.0/np.sqrt(1-X**2)
    return X*1.0/np.sqrt(1+X**2)
    #return ["%.3f" % z for z in V]

RPAlist = np.array([-30,-10, -5,-4,-3])

#for RPA in RPAlist:
#    ax2.axvline(RPA, linestyle="--", color='r')

#ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
ax1.yaxis.set_major_locator(MultipleLocator(0.1))

ax2.set_xticks(inv_tick_function(RPAlist))
ax2.set_xticklabels(np.abs(RPAlist))
ax2.set_xlim(ax1.get_xlim())

ax1.text(-0.942, 1.06, r"Ensemble " + expt, ha="right")
ax1.text(-0.9976, 0.53, r"Cross section discrepancy: $(\sigma^D - \sigma^M)^2/\sigma^D$", ha="left", fontsize=16.0)

#ax1.set_ylabel(r'$\sqrt{\lambda_p^{D \,2} + \lambda_p^{\overline{D}\, 2}}/\sqrt{\lambda_n^{D \,2} + \lambda_n^{\overline{D}\, 2}}$')
ax1.set_ylabel(r'$\lambda_n/\lambda_p$', fontsize=22.0)
ax1.set_xlabel(r'$f =  (\lambda_p^D \lambda_n^{D} + \lambda_p^{\overline{D}} \lambda_n^{\overline{D}})/ \sqrt{(\lambda_p^{D \,2} + \lambda_p^{\overline{D}\, 2})(\lambda_n^{D \,2} + \lambda_n^{\overline{D}\, 2})}$', fontsize=20)

ax2.set_xlabel(r'$|\lambda_n^D/\lambda_n^{\overline{D}}|$', fontsize=20.0)


#pl.figure()

#for i in range(32):
#    print i+1, "; l_n/l_p = ", getR(i+1)
#    print getSigvals(i+1)

#print pA
#sigs = getSigvals(pA)
#print sigs
#pl.hist(sigs)

#pl.tight_layout()

pl.savefig("../plots/Significance-Estimate.pdf", bbox_inches="tight")
pl.show()