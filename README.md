# TransversionsInRegElements

## Running the jaspar transistion / transversion analysis

This analysis investigates whether transversions are expected to have a bigger impact on 
TF binding affinity than transisions. The analysis relies on the JASPAR nonredundant motif set
to predift TF binding affinity, and a linear regression analysis to investigate the effects
of potential mutations while correcting for each motif and the position of the mutation
within each motif.

Performing the analysis is split across two scrips. The first script downloads the nonredundant jaspar
database and evaluate the predicted effect on TF affinity of every possible mutation at every 
possible position in every motif. This analysis is performed by the Effects_on_TFBS_motif_scores.py
script. 

The second step is the linear regression analyses to ask whether transversions are associated with
a bigger effect on TF binding affinity than transitions. That analysis is done in the jaspar_regression_analysis.R
script. The output is a log file of overall summary statistics for the data, regression statistics, plots
of the regression results, and an R image that canbe loaded for further invstigation. The R code requires that 
the reshape2 and car packages are installed.

Both steps of the analysis are combined in a simple shell script, jaspar_ts_tv_analysis.sh
