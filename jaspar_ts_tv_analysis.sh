#!/bin/bash

echo "Beginning jaspar Ts:Tv analysis"
python Effects_on_TFBS_motif_scores.py > jaspar_Ts_Tv.tab
echo "Calculated predicted effects of all possible mutations"
echo "Beginning regression analysis"
Rscript jaspar_regression_analysis.R
