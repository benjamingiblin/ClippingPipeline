#!/bin/bash

# 22/01/2021, B. Giblin, Postdoc, Edinburgh
# This code reads in two clipping pipeling param files, potentially corresponding to different
# redshift bins, then cycles through the simulation lines of sight,
# running TreeCorr to calculate the (cross-)correlation between them.
# It calculates unclippedXunclipped, clippedXclipped, and also unclippedXclipped correlations.

# This code relies on the clipping pipeline having already been ran, with the unclipped & clipped
# shear catalogues having been copied to the appropriate subdirectories on cuillin.

# This code is can be launched on cuillin with the script: Launch_Calc_CrossCorr_Cycle_z_And_LOS.sh
# You can also execute: Run_Launch_Calc_CrossCorr_Cycle_z_And_LOS.sh
# which cycles through the redshift bin combinations, launching the launch script.

# the parameter files specifying the redshift bin, clipping info etc.
paramfile1=$1
paramfile2=$2

# the lines of sight to cycle through
los_start=$3
los_end=$4

for los in `seq $los_start $los_end`; do
    python Tree_Correlation_Function/TreeCorr_CorrFun_CrossCorr.py Sims_Run $paramfile1 $los $los DUMMY Sims_Run $paramfile2 $los $los
done
