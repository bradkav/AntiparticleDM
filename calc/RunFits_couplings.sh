#!/bin/bash

#Experimental ensemble to use
EXPT=A

#DM mass in GeV
MASS=50

#Output directory for the results
OUTPUT_DIR=../results/new

mkdir $OUTPUT_DIR

#Loop over 1024 grid points in (f, l_n/l_p)
for i in 1..1024
do
    OUTFILE=$(OUTPUT_DIR)/Results_p$i.txt
    python CalcDisc-vs-Params.py $EXPT $MASS $i $OUTFILE
done