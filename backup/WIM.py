from DMUtils import *
from Experiment import *
import sys
#import matplotlib.pyplot as pl
from scipy.stats import chi2, norm
import CalcParamPoint as CPP
from CalcLikelihood import *

ensemble = sys.argv[1]
m0 = float(sys.argv[2])
index = int(sys.argv[3])

if (len(sys.argv) > 4):
    output_folder = sys.argv[4]
else:
    output_folder = "results/"


print " Loading experiments for ensemble", ensemble, "..."

if (ensemble == "A"):
    exptlist = ["Xenon2", "Argon", "Silicon"]
elif (ensemble == "B"):
    exptlist = ["Xenon2", "Argon", "Germanium"]
elif (ensemble == "C"):
    exptlist = ["Xenon2", "Argon", "CaWO4"]
elif (ensemble == "D"):
    exptlist = ["Xenon2", "Argon", "Germanium_half","CaWO4_half"]


#exptlist = ["Xenon2", "Argon", "Silicon"]
N_expt = len(exptlist)
expts = [ Experiment(exptlist[i] + ".txt") for i in range(N_expt)]

#print " Generating events..."

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
print " Couplings:", l
for i in range(N_expt):
    print " Ne(" + exptlist[i] + "): ", expts[i].CalcNevents(m0,l), "; sig_p =", expts[i].sig_eff(m0, l)
print " "
print " Calculating likelihoods..."

Nmvals = 20
mlist = np.logspace(np.log10(m0/2.0), np.log10(m0*2.0), Nmvals)

if (m0 > 200):
    mlist = np.logspace(np.log10(m0/10.0), np.log10(m0*10.0), Nmvals)

#print mlist
likelist_maj = np.zeros((Nmvals))
likelist_dir = np.zeros((Nmvals))

#Nlist = np.array((4, 25, 50, 100))

#100 is the correct number, but 50 is quick!

Nsamps = 50

sigvals = np.zeros(Nsamps)

#If we exceed Argon limits, just output zeroes...
if (expts[1].sig_eff(m0, l) > 1e-42):
    np.savetxt(output_folder + "Results_p" + str(index)+".txt",sigvals, header="Ensemble "+ ensemble + ", m = " + str(m0) + ", lambda = "+str(l) + "; ARGON LIMIT EXCEEDED")
    sys.exit()

for k in range(Nsamps):
    for expt in expts:
        expt.GenerateEvents(m0, l)
    for i, mi in enumerate(mlist):
        #print i+1, mi
        #for j in range(len(Nlist)):
        likelist_maj[i],likelist_dir[i] = CalcLike_grid(mi, expts, 100, refine=True)
    sig = CalcSignificance(np.nanmax(likelist_maj), np.nanmax(likelist_dir))
    sigvals[k] = sig
    print " Sample", str(k+1), " - Discrimination significance:", sig, "sigma"

print np.median(sigvals)


np.savetxt(output_folder + "Results_p" + str(index)+".txt",sigvals, \
    header="Ensemble "+ ensemble + ", m = " + str(m0) + ", lambda = "+str(l))
