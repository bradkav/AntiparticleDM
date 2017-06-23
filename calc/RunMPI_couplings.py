#!/usr/bin/env python
from mpi4py import MPI
from subprocess import call
import sys

#print sys.argv

expt = sys.argv[1]
mass = sys.argv[2]
arr_index = int(sys.argv[3])
resdir = sys.argv[4]

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

index = nprocs*(arr_index-1) + rank + 1

myDir = "/home/kavanagh/Antiparticle/WIMpy/"
cmd = "cd "+myDir+" ; python CalcDiscrimination.py "
cmd += expt + " " + mass + " " + str(index) + " " + resdir 
#print "cmd"
sts = call(cmd,shell=True)
comm.Barrier()
