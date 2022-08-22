#!/bin/bash

#SBATCH -n 2
#SBATCH -p all
#SBATCH -t 7-00:00
#SBATCH --job-name=CosmoSLICS_NoiseReal
#SBATCH --requeue
#SBATCH --mail-type=ALL
#SBATCH --constraint="3TBdatadisk"
#SBATCH --mem=15000

###### NB: Normally -mem=150000

###################### mem is in MB, based on watching memory usage of get_clipped_shear.py (most mem intensive) on eday.

cosmol=$1
los_start=$2
los_end=$3

./Master_CorrFun_ByParts.sh Sims_Run /home/bengib/Clipping_SimsLvW/param_files/100Sqdeg_SN0.28_NoMask_LSSTGpAM_X3sigma_SS17_zLSST_ZBcut0.6-1.4_ThBins9 $los_start $los_end

#### The 10arcsec res + shape noise SLICS used for Chris' project:
#### 100Sqdeg_SN0.28_NoMask_LSSTGpAM_X3sigma_SS17_zLSST_ZBcut0.6-1.4_ThBins9

#### The 60 arcsec noise-free SLICS currently being used in my PDF project:
#### 100Sqdeg_NF_NoMask_LSSTGpAM_X3sigma_SS1.41_zLSST_ZBcut0.6-1.4_ThBins9_MRres60arcsec

###########################SBATCH --mail-user=bengib@roe.ac.uk 
###########SBATCH -o log_Cosmol${i}.out 
