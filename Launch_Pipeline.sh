#!/bin/bash

#SBATCH -n 2
#SBATCH -p all
#SBATCH -t 7-00:00
#SBATCH --job-name=CosmoSLICS_NoiseReal
#SBATCH --requeue
#SBATCH --mail-type=ALL
#SBATCH --constraint=datadisk
#SBATCH --mem=15000

###### NB: Normally -mem=150000

###################### mem is in MB, based on watching memory usage of get_clipped_shear.py (most mem intensive) on eday.

./Master_CorrFun_ByParts.sh Sims_Run /home/sacreman/Clipping_Pipeline/param_files/100Sqdeg_SN0.27_NoMask_KiDS1000GpAM_X3sigma_SS9.33_zKiDS1000_ZBcut0.1-0.3_ThBins9_Cosmolfid_MRres60arcsec 1 50


###########SBATCH -o log_Cosmol${i}.out 
