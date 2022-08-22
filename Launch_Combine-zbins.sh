#!/bin/bash

#SBATCH -n 2
#SBATCH -p all
#SBATCH -t 7-00:00
#SBATCH --job-name=CosmoSLICS_NoiseReal
#SBATCH --requeue
#SBATCH --mail-type=ALL
#SBATCH --constraint="3TBdatadisk"
#SBATCH --mem=25000

paramfile1=$1
los_start=$2
los_end=$3
paramfile2=$4

./Master_CorrFun_ByParts.sh Sims_Run $paramfile1 $los_start $los_end $paramfile2


###########################SBATCH --mail-user=bengib@roe.ac.uk 
###########SBATCH -o log_Cosmol${i}.out 
