import sys
import CalcParamPoint as CPP
from CalcLikelihood import *
from WIMpy.Experiment import Experiment

print " "
print "**********************************"
print "*      CalcDiscrimination.py     *"
print "**********************************"
print " "

print " Code for calculating significance of discriminating between Dirac and Majorana DM..."
print " "

#Get inputs from command line
ensemble = sys.argv[1] # A, B, C or D
m0 = float(sys.argv[2]) # 10.0 - 10000.0
index = int(sys.argv[3]) # 1-1024

save_to_file = False

if (len(sys.argv) > 4):
    output_folder = sys.argv[4] #Absolute or relative path to output folder
    save_to_file = True

#Number of mock data sets to generate
Nsamps = 100

#Select the correct experimental ensemble
print " Loading experiments for ensemble", ensemble, "..."
if (ensemble == "A"):
    exptlist = ["Xenon", "Argon", "Silicon"]
elif (ensemble == "B"):
    exptlist = ["Xenon", "Argon", "Germanium"]
elif (ensemble == "C"):
    exptlist = ["Xenon", "Argon", "CaWO4"]
elif (ensemble == "D"):
    exptlist = ["Xenon", "Argon", "Germanium_half","CaWO4_half"]

exptdir = "DDexpt/"

#Initialise experiments
N_expt = len(exptlist)
expts = [ Experiment(exptdir + exptlist[i] + ".txt") for i in range(N_expt)]

#Generate couplings from index:
l = CPP.GeneratePoint_ind(index)

#This is where we calculate the normalisation of the couplings
#Aim for the same number of events in Xenon as:
# sigma = 1e-46 and m = 50
sig0 = expts[0].sig_eff(50, l)
l *= np.sqrt(1e-46/sig0)
Ntarget = expts[0].CalcNevents(50, l)
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


#List of masses to calculate for
Nmvals = 20
mlist = np.logspace(np.log10(m0/2.0), np.log10(m0*2.0), Nmvals)

if (m0 > 200):
    mlist = np.logspace(np.log10(m0/10.0), np.log10(m0*10.0), Nmvals)

likelist_maj = np.zeros((Nmvals))
likelist_dir = np.zeros((Nmvals))


print " Calculating likelihoods..."

#Discrimination significance
sigvals = np.zeros(Nsamps)

#If we exceed Argon limits, just output zeroes...
if (expts[1].sig_eff(m0, l) > 1e-42):
    print " Argon limit exceeded..."
    if (save_to_file):
        np.savetxt(output_folder + "Results_p" + str(index)+".txt",sigvals, header="Ensemble "+ ensemble + ", m = " + str(m0) + ", lambda = "+str(l) + "; ARGON LIMIT EXCEEDED")
    sys.exit()

#Generate Nsamps mock data sets
for k in range(Nsamps):
    #Generate sample of events
    for expt in expts:
        expt.GenerateEvents(m0, l)
    #Calculate likelihood on a grid for range of mass values
    for i, mi in enumerate(mlist):
        #Majorana hypothesis
        likelist_maj[i] = CalcLike_grid(mi, expts, 200, maj=True, refine=True)
        #Dirac hypothesis
        likelist_dir[i] = CalcLike_grid(mi, expts, 50, maj=False, refine=True)
    sig = CalcSignificance(np.nanmax(likelist_maj), np.nanmax(likelist_dir))
    sigvals[k] = sig
    print " Sample", str(k+1), " - Discrimination significance:", sig, "sigma"

print " Median significance:", np.median(sigvals)

#Output to file
if (save_to_file):
    np.savetxt(output_folder + "Results_p" + str(index)+".txt",sigvals, \
        header="Ensemble "+ ensemble + ", m = " + str(m0) + ", lambda = "+str(l))
