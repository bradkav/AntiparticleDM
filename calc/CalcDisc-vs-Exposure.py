"""
CalcDisc-vs-Exposure.py

This file is a wrapper for CalcDiscrimination.py which
allows you to specify the exposure of the experiments
involved (keeping the Xenon and Argon exposures fixed),
and calculates the discrimination significance in that case.

BJK - 23/06/2017
"""

import sys
import numpy as np
import CalcParamPoint as CPP
import CalcDiscrimination as CD

#Get inputs from command line
ensemble = sys.argv[1] # ensemble: A, B, C or D
m0 = float(sys.argv[2]) # DM mass [GeV]: 10.0 - 10000.0
r_np = float(sys.argv[3]) # Value of lambda_n/lambda_p
expindex = int(sys.argv[4])-1 # Index specifying the exposure (see below): 1-32

if (len(sys.argv) > 5):
    outfile = sys.argv[5] #Absolute or relative path to output file
else:
    outfile = None

#Fix the value of f:
f = -0.995

#Calculate exposure in kg-years  
exprange = np.round(np.logspace(2, 5 , 32))                                   
exp = int(exprange[expindex])

#Calculate discrimination significance
CD.CalcDiscrim(ensemble, m0, f, r_np, outfile, exposure=exp)

#Output to file
#if (save_to_file):
#    np.savetxt(output_folder + "Results_" + str(exp) + "kg.txt",sigvals, \
#        header="Ensemble "+ ensemble + ", m = " + str(m0) + ", c_n/c_p = " + R0 + ", lambda = "+str(l)\
#        + ", exposure = " + str(exp) + "kg years (before eff.)")
