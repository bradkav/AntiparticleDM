#!/usr/bin/env python
from mpi4py import MPI
from subprocess import call
import sys

#print sys.argv

expt = sys.argv[1]
mass = sys.argv[2]
Rval = sys.argv[3]
arr_index = int(sys.argv[4])
resdir = sys.argv[5]

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

index = nprocs*(arr_index-1) + rank + 1

myDir = "/home/kavanagh/Antiparticle/WIMpy/"
cmd = "cd "+myDir+" ; python CalcDisc-vs-Exposure.py "
cmd += expt + " " + mass + " " + Rval + " " + str(index) + " " + resdir 
#print "cmd"
sts = call(cmd,shell=True)
comm.Barrier()
