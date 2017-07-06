#!/bin/bash

for EXPT in A B C D
do 
	for MASS in 25 50 300 1000
	do
		python ../analysis/PlotContours.py $EXPT $MASS
	done
done