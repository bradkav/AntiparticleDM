"""
CalcDisc-vs-Couplings.py

This file is a wrapper for CalcDiscrimination.py which
allows you to specify an index 1-1024 (which denotes
a particular point in parameter space) and then
calculates the discrimination significance at that point.

BJK - 23/06/2017
"""

import sys
import CalcParamPoint as CPP
import CalcDiscrimination as CD

#Get inputs from command line
ensemble = sys.argv[1] # ensemble: A, B, C or D
m0 = float(sys.argv[2]) # DM mass [GeV]: 10.0 - 10000.0
index = int(sys.argv[3]) # Index of param point: 1-1024

print " Calculating for index: ", index

if (len(sys.argv) > 4):
    outfile = sys.argv[4] #Absolute or relative path to output file
else:
    outfile = None
    
#Calculate parameter values from index
fval = CPP.getf(index)
Rval = CPP.getR(index)
    
#Calculate discrimination significance
CD.CalcDiscrim(ensemble, m0, fval, Rval, outfile)
    
##Output to file
#if (save_to_file):
#    np.savetxt(output_folder + "Results_p" + str(index)+ ".txt",sigvals, \
#        header="Ensemble "+ ensemble + ", m = " + str(m0) + ", lambda = "+str(l))
