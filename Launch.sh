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

#cosmol=$1
SN=$1
ZBcut=$2
los_start=$3
los_end=$4

./Master_CorrFun_ByParts.sh Sims_Run /home/bengib/Clipping_SimsLvW/param_files/100Sqdeg_SN${SN}_NoMask_KiDS1000GpAM_X3sigma_SS1_zKiDS1000_ZBcut${ZBcut}_ThBins9_MRres60arcsec $los_start $los_end


#### The 10arcsec res + shape noise SLICS used for Chris' project:
#### 100Sqdeg_SN0.28_NoMask_LSSTGpAM_X3sigma_SS17_zLSST_ZBcut0.6-1.4_ThBins9

#### The 60 arcsec noisy K1000-like SLICS/cosmoSLICS currently used in my PDF project
#### 100Sqdeg_SN0.273_NoMask_KiDS1000GpAM_X3sigma_SS1_zKiDS1000_ZBcut0.5-0.7_ThBins9_Cosmolfid_MRres60arcsec

#### The 140.625 arcsec noise-free (smoothing scale 2.3 arcmin) current being used for Rhys' CNN
#### param_files/100Sqdeg_NF_NoMask_3.32GpAM_X3sigma_SS1_zKiDS_ZBcut0.9-1.2_ThBins9_MRres140.625arcsec

###########################SBATCH --mail-user=bengib@roe.ac.uk 
###########SBATCH -o log_Cosmol${i}.out 
