from DMUtils import *
from Experiment import *
import sys
#import matplotlib.pyplot as pl
from scipy.stats import chi2, norm
import CalcParamPoint as CPP

ensemble = sys.argv[1]
m0 = float(sys.argv[2])
index = int(sys.argv[3])


if (len(sys.argv) > 4):
    output_folder = sys.argv[4]
else:
    output_folder = "results/"

target_sigma = 1e-46
if (m0 < 30):
    target_sigma = 2e-46

if (m0 > 90):
    target_sigma = 1e-45

#Recalculate to give the correct numbers!
l = CPP.GeneratePoint_ind(index)

#l = CPP.GeneratePoint(1/0.8, -0.992)
sig0 = (1.973e-14*1.973e-14)*4.0*(reduced_m(1.0, m0))**2.0/np.pi
sig = sig0*0.5*((l[0]*54.0 + l[1]*77.0)**2.0 + (l[2]*54.0 + l[3]*77.0)**2.0)
sig_p = sig/(54+77)**2
l *= np.sqrt(target_sigma/sig_p)

#l = [  1.01840302e-07 , -7.26022796e-08 ,  0.00000000e+00  , 0.00000000e+00]
#print l
#Should be 312.544 events!

#Sampling paramaters
loglsq_min = -11
loglsq_max = -7

Nvals = 50

#----Functions----


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

for i in range(N_expt):
    print " Ne(" + exptlist[i] + "): ", expts[i].CalcNevents(m0,l), "; sig_p =", expts[i].sig_eff(m0, l)

sigvals = np.zeros((32,32))
fvals = 
for i in range(32*32):
    l = CPP.GeneratePoint_ind(i+1)
    

