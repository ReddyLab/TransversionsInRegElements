#!/bin/bash

python Effects_on_TFBS_motif_scores.py > jaspar_Ts_Tv.tab
Rscript jaspar_regression_analysis.R
