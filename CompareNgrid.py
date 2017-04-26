from DMUtils import *
from Experiment import *
import sys
import matplotlib.pyplot as pl
from scipy.stats import chi2, norm
import emcee
import CalcParamPoint as CPP
from CalcLikelihood import *


print " "
print "***********************"
print "*   CompareNgrid.py   *"
print "***********************"
print " "

print " Code for comparing sampling grids for calculating the likelihood..."



#Read parameters from command line
m0 = float(sys.argv[1])
index = int(sys.argv[2])


#----Functions----

print " Loading experiments..."
exptlist = ["Xenon2", "Argon", "Silicon"]
N_expt = len(exptlist)
expts = [ Experiment(exptlist[i] + ".txt") for i in range(N_expt)]


#Generate couplings from index:                                               
l = CPP.GeneratePoint_ind(index)

#This is where we calculate the normalisation of the couplings                
#Aim for the same number of events in Xenon as:                               
# sigma = 1e-46 and m = 50                                                    
sig0 = expts[0].sig_eff(50, l)
l *= np.sqrt(1e-46/sig0)
Ntarget = expts[0].CalcNevents(50, l)
print " Target number of events in Xenon:", Ntarget

#Rescale to give correct event numbers for mass m0                            
l *= np.sqrt(Ntarget/expts[0].CalcNevents(m0, l))


#l is in the format [lpD, lnD, lpDb, lnDb]
print " "
print " DM mass [GeV]:", m0
print " lambda [GeV^-2]:", l
print " f =", CPP.getf(index)
print " c_n/c_p =", CPP.Calc_Rnp(l)
print " "
for i in range(N_expt):
    print " Ne(" + exptlist[i] + "): ", expts[i].CalcNevents(m0,l), "; sig_p =", expts[i].sig_eff(m0, l)
print " "


print " Generating events..."
for i in range(N_expt):
    expts[i].GenerateEvents(m0, l)
    print " Nobs(" + exptlist[i] + "): ",len(expts[i].events)
print " "

print " Calculating likelihoods..."
Nmvals = 20

m_min = m0/2.0
m_max = m0*2.0
if (m0 > 200):
    m_min = m0/10.0
    m_max = m0*10.0

mlist = np.logspace(np.log10(m_min), np.log10(m_max), Nmvals)



likelist = np.zeros((Nmvals, 5))

#Different numbers of grid points to try
Ngrid = [50, 100, 200]
refine = [True, True, True]
#In the last case, we use 100 points but refine the grid

likelist_maj = np.zeros((Nmvals, 3))
likelist_dir = np.zeros((Nmvals, 3))

for i, mi in enumerate(mlist):
    print "   ",i+1, "of", Nmvals,": m_x =", mi, "GeV"
    for j in range(1):
        #print j
        likelist_maj[i,j] = CalcLike_grid(mi, expts, 4*Ngrid[j], maj = True,refine=refine[j])
        likelist_dir[i,j] = CalcLike_grid(mi, expts, Ngrid[j], maj = False,refine=refine[j])
        if (likelist_dir[i,j] < likelist_maj[i,j]):
            likelist_dir[i,j] = likelist_maj[i,j]
        
for j in range(3):
    sig = CalcSignificance(np.nanmax(likelist_maj[:,j]), np.nanmax(likelist_dir[:,j]))
    str_extra = ""
    if (refine[j]):
        str_extra = ", refined"
    print " Discrimination significance (N_grid = "+str(Ngrid[j])+str_extra+"):", sig, "sigma"


lines = [":", "--", "-"]
labels = ["N = 50 (refined)",  "N = 100 (refined)", "N = 200 (refined)"]

pl.figure()
L0_global = np.nanmax(likelist_dir)
for j in range(3):
    L0 = np.nanmax(likelist_dir[:,j])
    #print L0
    pl.semilogx(mlist, -2*(likelist_maj[:,j]-L0_global), 'b',\
        linestyle=lines[j], linewidth=1.5)
    pl.semilogx(mlist, -2*(likelist_dir[:,j]-L0_global), 'r',\
        linestyle=lines[j], linewidth=1.5)

#Add dummy lines for labels
pl.semilogx(1e-30, 1e-30, 'r-',label=r"Dirac", linewidth=1.5)
pl.semilogx(1e-30, 1e-30, 'b-',label=r"Majorana", linewidth=1.5)
for j in range(3):
    pl.semilogx(1e-30, 1e-30, 'k',\
         linestyle=lines[j],label=labels[j], linewidth=1.5)


pl.legend(loc="best", frameon=False)
pl.ylim(-1, 50)
pl.xlim(m_min, m_max)
pl.axvline(m0, linestyle='--', color='k')
pl.axhline(0, linestyle='--', color='k')
pl.xlabel(r"$m_\chi [GeV]$")
pl.ylabel(r"$-2 \Delta \mathrm{log}\mathcal{L}$")


#pl.savefig("plots/GridComparison_m=" + str(m0) +".pdf")
pl.show()

