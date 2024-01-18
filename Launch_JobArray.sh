#!/bin/bash

#SBATCH --array=1-100
#SBATCH -n 2
#SBATCH -p all
#SBATCH -t 7-00:00
#SBATCH --job-name=CosmoSLICS_PDFs_0-9
#SBATCH --requeue
#SBATCH --mail-type=ALL
#SBATCH --constraint="3TBdatadisk"
#SBATCH --mem=25000

# Specify the path to the config file
config=config_cross2_RAND.txt 

# Extract the paramfile for the current $SLURM_ARRAY_TASK_ID
paramfile=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
# Extract the los_start & end
los_start=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
los_end=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)
# 2nd paramfile
paramfile2=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $config)

# TEST TO MAKE SURE JOB ARRAY WORKS:
# Print to a file a message that includes the current $SLURM_ARRAY_TASK_ID,
# the paramfile, los_start and los_end 
#echo "This is array task ${SLURM_ARRAY_TASK_ID}, the paramfile is ${paramfile} and los_start,end is ${los_start},${los_end}." >> output${SLURM_ARRAY_TASK_ID}.txt

#python Correlation_Function/plot_CorrFun.py Sims_Run $paramfile $los_start $los_end
#./Master_CorrFun_ByParts.sh Sims_Run $paramfile $los_start $los_end
./Master_CorrFun_ByParts.sh Sims_Run $paramfile $los_start $los_end $paramfile2


###########################SBATCH --mail-user=bengib@roe.ac.uk 
###########SBATCH -o log_Cosmol${i}.out 
