#!/bin/bash

#Experimental ensemble to use
EXPT=A

#DM mass in GeV
MASS=50

#Coupling values (f = -0.995 is fixed in CalcDisc-vs-Exposure.py)
#Ratio of neutron to proton couplings
RNP=0.75

#Output directory for the results
OUTDIR=../results/new_exposure

mkdir -p $OUTDIR

#Loop over 20 different values of the exposure
trap "exit" INT
for i in {1..20}
do
    #Choose whatever naming scheme you want for the output files
    OUTFILE=$OUTDIR/Results_r${RNP}_exp$i.txt

    #Calculate significance
    python CalcDisc-vs-Exposure.py $EXPT $MASS $RNP $i $OUTFILE
done