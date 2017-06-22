import sys
import CalcParamPoint as CPP
from CalcLikelihood import *
from WIMpy.Experiment import Experiment
from scipy.interpolate import interp1d
import matplotlib.pyplot as pl

print " "
print "**********************************"
print "*      CalcDisc-vs-Exposure.py     *"
print "**********************************"
print " "

print " Code for calculating significance of discriminating between Dirac and Majorana DM..."
print " "

#Get inputs from command line
ensemble = sys.argv[1] # A, B, C or D
m0 = float(sys.argv[2]) # 10.0 - 10000.0
R0 = float(sys.argv[3]) # Value of c_n/c_p
expindex = int(sys.argv[4])-1 # 1-20

#Calculate exposure       
exprange = np.round(np.logspace(2, 5 , 20))                                   
exp = int(exprange[expindex])

save_to_file = False

if (len(sys.argv) > 5):
    output_folder = sys.argv[5] #Absolute or relative path to output folder
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

#Set exposure...                                                              

#Correction for number of experiments                                         
corr = 1.0/(N_expt - 2.0)

for i in range(2, N_expt):
    expts[i].exposure = corr*exp*365*0.7

print " exposure = "+str(expts[-1].exposure/365.0) + " kg yr"


#Generate couplings from index:
l = CPP.GeneratePoint(R0, -0.995)

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
print " f = - 0.995"
print " c_n/c_p =", CPP.Calc_Rnp(l)
print " "
for i in range(N_expt):
    expts[i].TabulateAll(m0)
    print " Ne(" + exptlist[i] + "): ", expts[i].CalcNevents(m0,l), "; sig_p =", expts[i].sig_eff(m0, l)
print " "


#If we exceed Argon limits, just output zeroes...
if (expts[1].sig_eff(m0, l) > 1e-43):
    np.savetxt(output_folder + "Results_p" + str(index)+".txt",sigvals, header="Ensemble "+ ensemble + ", m = " + str(m0) + ", lambda = "+str(l) + "; ARGON LIMIT EXCEEDED")
    sys.exit()


#List of masses to calculate for
Nmvals = 25

delta_m = np.log10(m0-10.0)

m_min = m0/delta_m
m_max = m0*delta_m
if (m_min < 20.0):
    m_min = 20.0


if (m0 >= 500):
    m_min = m0/10.0
    m_max = m0*10.0
    
mlist = np.logspace(np.log10(m_min), np.log10(m_max), Nmvals)

likelist_maj = np.zeros((Nmvals))
likelist_dir = np.zeros((Nmvals))


print " Calculating likelihoods..."

#Discrimination significance
sigvals = np.zeros(Nsamps)

#Generate Nsamps mock data sets
for k in range(Nsamps):
    #Generate sample of events
    for expt in expts:
        expt.GenerateEvents(m0, l)
    #Calculate likelihood on a grid for range of mass values
    for i, mi in enumerate(mlist):
        #Majorana hypothesis
        #print mi
        likelist_maj[i] = CalcLike_grid(mi, expts, 200, maj=True, refine=True)
        #Dirac hypothesis
        likelist_dir[i] = CalcLike_grid(mi, expts, 50, maj=False, refine=True)
    sig = CalcSignificance(np.nanmax(likelist_maj), np.nanmax(likelist_dir))
    sigvals[k] = sig
    print " Sample", str(k+1), " - Discrimination significance:", sig, "sigma"
    #L0_global = np.nanmax(likelist_dir)
    #pl.figure()
    #pl.plot(mlist, -2*(likelist_maj - L0_global), 'b--')
    #pl.plot(mlist, -2*(likelist_dir - L0_global), 'r--')
    #pl.ylim(-1, 60)
    #pl.show()

print " Median significance:", np.median(sigvals)

#Output to file
if (save_to_file):
    np.savetxt(output_folder + "Results_" + str(exp) + "kg.txt",sigvals, \
        header="Ensemble "+ ensemble + ", m = " + str(m0) + ", c_n/c_p = " + R0 + ", lambda = "+str(l)\
        + ", exposure = " + str(exp) + "kg years (before eff.)")
