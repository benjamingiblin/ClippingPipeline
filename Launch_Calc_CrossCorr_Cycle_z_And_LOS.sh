#!/bin/bash

#SBATCH -n 2
#SBATCH -p all
#SBATCH -t 7-00:00
#SBATCH --job-name=CosmoSLICS_NoiseReal
#SBATCH --requeue
#SBATCH --mail-type=ALL
#SBATCH --constraint="3TBdatadisk"
#SBATCH --mem=15000

# the parameter files specifying the redshift bin, clipping info etc.
paramfile1=$1
paramfile2=$2

# the lines of sight to cycle through
los_start=$3
los_end=$4

R_start=$5
R_end=$6

#./Calc_CrossCorr_Cycle_z_And_LOS.sh $paramfile1 $paramfile2 $los_start $los_end $R_start $R_end
python Correlation_Function/plot_CorrFun.py Sims_Run $paramfile1 $los_start $los_end $paramfile2


###########SBATCH -o log_Cosmol${i}.out 
