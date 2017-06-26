#!/usr/bin/env python
from mpi4py import MPI
from subprocess import call
import sys

"""
This script will create 'nprocs' parallel instances of 
CalcDisc-vs-Couplings.py and calculate the significance
for a different coupling index on each instance.

For a given mass and ensemble, there are 1024 indices, so
this code is meant to be run in batches. The number of the
batch is specified in arr_index (starting from arr_index=1).
The code will then automatically calculate the appropriate
coupling indices based on the batch number and the number of
processes in the current batch.

See CalcDisc-vs-Couplings.py for information on what the other
command line parameters should be.

To execute with 16 MPI processes, you would have something
like the following in your submission script:

mpirun -np 16 python RunMPI_couplings.py

"""

#Get inputs from command line
ensemble = sys.argv[1]
mass = sys.argv[2]
arr_index = int(sys.argv[3]) #Array index of the current batch
output_dir = sys.argv[4] #output directory for the results

comm = MPI.COMM_WORLD
#Get total number of MPI processes
nprocs = comm.Get_size()
#Get rank of current process
rank = comm.Get_rank()

#Coupling index
index = nprocs*(arr_index-1) + rank + 1
#Output file labelled by coupling index
outfile = output_dir + "Results_p" + str(index) + ".txt"

#Directory where the calc files are located
myDir = "/home/kavanagh/AntiparticleDM/calc/"
cmd = "cd "+myDir+" ; python CalcDisc-vs-Couplings.py "
cmd += ensemble + " " + mass + " " + str(index) + " " + outfile 

sts = call(cmd,shell=True)
comm.Barrier()
