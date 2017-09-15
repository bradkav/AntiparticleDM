#!/bin/bash

python ../analysis/PlotContours_row.py A	#Fig. 1
python ../analysis/PlotContours_row.py B	#Fig. 2
python ../analysis/PlotContours_row.py C	#Fig. 3
python ../analysis/PlotContours_row.py D	#Fig. 4

python ../analysis/PlotExposure.py 0.75	#Fig. 5 (left)
python ../analysis/PlotExposure.py 0.8	#Fig. 5 (right)

python ../analysis/PlotFundamentalCouplings.py 	#Fig. 6

python ../analysis/PlotContours_analytic.py 	#Fig. 7

python ../analysis/PlotExposure_comparison.py 0.75	#Comparison of single and multiple isotope calculations 