#!/usr/bin/env python
from mpi4py import MPI
from subprocess import call
import sys

"""
This script will create 'nprocs' parallel instances of 
CalcDisc-vs-Exposure.py and calculate the significance
for a different exposure on each instance.

For a given mass and ensemble, there are 32 exposures to
calculate for, so this code is meant to be run in batches. 
The number of the batch is specified in arr_index (starting
from arr_index=1). The code will then automatically calculate 
the appropriate index which labels the exposure to be used,
based on the batch number and the number of processes in 
the current batch.

See CalcDisc-vs-Exposure.py for information on what the other
command line parameters should be.
"""

#Get inputs from command line
expt = sys.argv[1] #First three arguments are as in CalcDisc-vs-Exposure.py
mass = sys.argv[2]
r_np = sys.argv[3]
arr_index = int(sys.argv[4]) #Array index of the current batch
output_dir = sys.argv[5] #output directory for the results

comm = MPI.COMM_WORLD
#Get total number of MPI processes
nprocs = comm.Get_size()
#Get rank of current process
rank = comm.Get_rank()

#Index labelling the exposure to use
expindex = nprocs*(arr_index-1) + rank + 1
#Output file labelled by r_np and exposure index
outfile = output_dir + "Results_r" + str(r_np) + "_exp" + str(expindex) + ".txt"

#Directory where the calc files are located
myDir = "/home/kavanagh/AntiparticleDM/calc/"
cmd = "cd "+myDir+" ; python CalcDisc-vs-Exposure.py "
cmd += expt + " " + mass + " " + r_np + " " + str(expindex) + " " + outfile 

sts = call(cmd,shell=True)
comm.Barrier()
