#!/bin/bash

#Experimental ensemble to use
EXPT=A

#DM mass in GeV
MASS=50

#Output directory for the results
OUTDIR=../results/new

mkdir -p $OUTDIR

#Loop over 1024 grid points in (f, l_n/l_p)
trap "exit" INT
for i in {1..1024}
do
    #Choose whatever naming scheme you want for the output files
    OUTFILE=$OUTDIR/Results_p$i.txt

    #Calculate significance
    python CalcDisc-vs-Couplings.py $EXPT $MASS $i $OUTFILE
done