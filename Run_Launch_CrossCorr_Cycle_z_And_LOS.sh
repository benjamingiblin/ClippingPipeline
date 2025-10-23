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


MRres="140.64arcs" # 60arcsec
Mask="Mosaic"
SS=2.816 #3.11, 9.33, 18.66  

# Redshift bins to cycle through, with potentially different SN levels
ZBcut=("0.1-0.3" "0.3-0.5" "0.5-0.7" "0.7-0.9" "0.9-1.2")

Survey="KiDS1000"
SN=(0.27 0.258 0.273 0.254 0.27)

dz_vals=("dz2" "dz3" "dz4" "dz5")

#Survey="LSST"
#SN=(0.28 0.28 0.28 0.28 0.28)

# The los to cycle through
los_start=74 #74 #74 #1
los_end=74 #292 #699 #50
R_start=1
R_end=18

missing_los=(135 140 449 595 596 597 598 601 610 614 735)

param_dir=/home/bengib/Clipping_Pipeline/param_files
# launch jobs for each cosmology and redshift bins:

count=1
#for i in `seq 0 24`; do         # Use this line for cosmoSLICS
for i in fid; do                 # this line for SLICS (and edit paramfile1/2 below).
    echo " ----------------------------------- On cosmol $i ----------------------------------- "
    
    for j in `seq 0 4`; do
	jp1=$((j+1))                 # use jp1 in line below if you want to avoid auto-bins
	for k in `seq $jp1 4`; do    

	    for dz1 in `seq 0 3`; do for dz2 in `seq 0 3`; do
					 if [ "$dz1" -eq "$dz2" ]; then
					     Sys_Tag1=${dz_vals[$dz1]}_; Sys_Tag2=$Sys_Tag1;
					 else
					     Sys_Tag1=${dz_vals[$dz1]}_; Sys_Tag2=${dz_vals[$dz2]}_;
					 fi  


					 
	    # !!! SET TO dz MOCKS !!!
	    paramfile1=$param_dir/${Sys_Tag1}100Sqdeg_SN${SN[$j]}_${Mask}_${Survey}GpAM_X3sigma_SS${SS}_z${Survey}_ZBcut${ZBcut[$j]}_ThBins9_Cosmol${i}_MRres${MRres}ec
	    paramfile2=$param_dir/${Sys_Tag2}100Sqdeg_SN${SN[$k]}_${Mask}_${Survey}GpAM_X3sigma_SS${SS}_z${Survey}_ZBcut${ZBcut[$k]}_ThBins9_Cosmol${i}_MRres${MRres}ec


	    echo "Running zbin $((j+1)) X zbin $((k+1))"
	    #ls $paramfile1
	    #ls $paramfile2
	    #echo "$count $paramfile1 $paramfile2 $los_start $los_end $R_start $R_end" >> config_corrfun.txt
	    sbatch Launch_Calc_CrossCorr_Cycle_z_And_LOS.sh $paramfile1 $paramfile2 $los_start $los_end $R_start $R_end
	    sleep 3s
	    #sbatch Launch_Combine-zbins.sh $i ${SN[$j]} ${ZBcut[$j]} ${SN[$k]} ${ZBcut[$k]} $SS $los_start $los_end
	    #python Correlation_Function/plot_CorrFun.py Sims_Run $paramfile1 $los_start $los_end $paramfile2
	    count=$((count+1))
	    
	    done; done
	    
	done
    done
    
done

echo $count

