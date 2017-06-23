#!/bin/bash

python PlotContours_row.py A	#Fig. 1
python PlotContours_row.py B	#Fig. 2
python PlotContours_row.py C	#Fig. 3
python PlotContours_row.py D	#Fig. 4

python PlotExposure.py 0.75	#Fig. 5 (left)
python PlotExposure.py 0.8	#Fig. 5 (right)

python PlotContours_analytic.py $Fig. 6