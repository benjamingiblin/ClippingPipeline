#!/bin/bash

#SBATCH -n 1
####SBATCH --cpus-per-task 8
#SBATCH -p all
#SBATCH -t 7-00:00
#SBATCH --job-name=CosmoSLICS_NoiseReal
#SBATCH --requeue
#SBATCH --mail-type=ALL
#SBATCH --mem=15000

###### NB: Normally -mem=150000

cosmology=$1
ia=$2
python Make_SNR-PDF_4cosmoSLICS.py $cosmology 
#python Make_SNR-PDF_4cosmoSLICS_NonMosaic.py $cosmology $ia

