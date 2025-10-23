#!/bin/bash

#SBATCH -n 2
#SBATCH -p all
#SBATCH -t 7-00:00
#SBATCH --job-name=CosmoSLICS_NoiseReal
#SBATCH --requeue
#SBATCH --mail-type=ALL
#SBATCH --constraint=datadisk
#SBATCH --exclude=worker[082,096]
#SBATCH --mem=25000

paramfile=$1
los_start=$2
los_end=$3
paramfile2=$4

#python Tree_Correlation_Function/TreeCorr_CorrFun.py Sims_Run $paramfile $los_start $los_end
#python Correlation_Function/plot_CorrFun.py Sims_Run $paramfile $los_start $los_end $paramfile2

# NOTE THIS SAYS Sims_Run!!!
./Master_CorrFun_ByParts.sh Sims_Run $paramfile $los_start $los_end


###########################SBATCH --mail-user=bengib@roe.ac.uk 
###########SBATCH -o log_Cosmol${i}.out 
