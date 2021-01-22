#!/bin/bash

# 22/01/2021, B. Giblin, Postdoc, Edinburgh
# This code runs cycles through the redshift bin combinations, launching jobs
# that calculate the unclipped & clipped cross-correlations between the redshifts.

# File hierarchy:
# Run_Launch_CrossCorr_Cycle_z_And_LOS.sh
# ...executes:
# Launch_Calc_CrossCorr_Cycle_z_And_LOS.sh
# ...which launches:
# Calc_CrossCorr_Cycle_z_And_LOS.sh
# this last code cycles through the simulation LOS, running this TreeCorr code at each step:
# Tree_Correlation_Function/TreeCorr_CorrFun_CrossCorr.py



# Redshift bins to cycle through, with potentially different SN levels
ZBcut=("0.1-0.3" "0.3-0.5" "0.5-0.7" "0.7-0.9" "0.9-1.2")
SN=(0.27 0.258 0.273 0.254 0.27)

# The los to cycle through
los_start=1
los_end=50

SS=9.33 # the smoothing scale used in clipping

param_dir=/home/bengib/Clipping_SimsLvW/param_files
# launch jobs for each cosmology and redshift bins:

for i in `seq 0 24`; do         # Use this line for cosmoSLICS
#for i in fid; do                 # this line for SLICS (and edit paramfile1/2 below).
    echo " ----------------------------------- On cosmol $i ----------------------------------- "
    
    for j in `seq 0 4`; do
	for k in `seq 0 4`; do

	    paramfile1=$param_dir/100Sqdeg_SN${SN[$j]}_NoMask_KiDS1000GpAM_X3sigma_SS${SS}_zKiDS1000_ZBcut${ZBcut[$j]}_ThBins9_Cosmol${i}_MRres60arcsec
	    paramfile2=$param_dir/100Sqdeg_SN${SN[$k]}_NoMask_KiDS1000GpAM_X3sigma_SS${SS}_zKiDS1000_ZBcut${ZBcut[$k]}_ThBins9_Cosmol${i}_MRres60arcsec

	    if [ "$k" -ge "$j" ]; then
		# (if statement allows bin1-bin2 to be ran, but not bin2-bin1)
		echo "Running zbin $((j+1)) X zbin $((k+1))"
		#ls $paramfile1
		#ls $paramfile2
		sbatch Launch_Calc_CrossCorr_Cycle_z_And_LOS.sh $paramfile1 $paramfile2 $los_start $los_end
	    fi
	    
		
	done
    done
    
done


