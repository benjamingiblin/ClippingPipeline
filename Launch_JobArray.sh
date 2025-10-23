#!/bin/bash

#SBATCH --array=1-219
#SBATCH -n 2
#SBATCH -p all
#SBATCH -t 7-00:00
#SBATCH --job-name=KiDS-NOISE_NoTomo
#SBATCH --requeue
#SBATCH --mail-type=ALL
#SBATCH --constraint=datadisk
#SBATCH --exclude=worker[082,096] 
#SBATCH --mem=25000

# ! GIBLIN ! MAKE SURE YOU CHANGE THE ARRAY NUM AT TOP TO NUM-ROWS IN CONFIG

# Specify the path to the config file
config=config_NOISE_SS2.816.txt

# Extract the paramfile for the current $SLURM_ARRAY_TASK_ID
paramfile=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
# Extract the los_start & end
los_start=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
los_end=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)
# 2nd paramfile
paramfile2=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $config)

# JOB OPTIONS:
#python Correlation_Function/plot_CorrFun.py Sims_Run $paramfile $los_start $los_end $paramfile2

#echo $paramfile "..." $los_start "..." $los_end
#./Master_CorrFun_ByParts.sh Sims_Run $paramfile $los_start $los_end
./Master_CorrFun_ByParts.sh Sims_Run $paramfile $los_start $los_end $paramfile2 

# this is the one for Launch_Calc_CrossCorr_Cycle_z_And_LOS.sh
#paramfile=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
#paramfile2=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
#los_start=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)
#los_end=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $config)

#R_start=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $6}' $config)
#R_end=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $7}' $config)
#n_start=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $8}' $config)
#n_end=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $9}' $config)
#./Calc_CrossCorr_Cycle_z_And_LOS.sh $paramfile $paramfile2 $los_start $los_end $R_start $R_end $n_start $n_end



# TEST TO MAKE SURE JOB ARRAY WORKS:
# Print to a file a message that includes the current $SLURM_ARRAY_TASK_ID,
# the paramfile, los_start and los_end 
#echo "This is array task ${SLURM_ARRAY_TASK_ID}, the paramfile is ${paramfile} and los_start,end is ${los_start},${los_end}." >> output${SLURM_ARRAY_TASK_ID}.txt


