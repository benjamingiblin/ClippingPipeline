#!/bin/bash

#SBATCH -n 2
#SBATCH -p all
#SBATCH -t 7-00:00
#SBATCH --job-name=CosmoSLICS_NoiseReal
#SBATCH --requeue
#SBATCH --mail-type=ALL
#SBATCH --constraint="3TBdatadisk"
#SBATCH --mem=25000

###cosmol=$1
###SN=$2
###ZBcut=$3
###SS=$4
###los_start=$5
###los_end=$6
##param_dir=/home/bengib/Clipping_SimsLvW/param_files
##paramfile=$param_dir/100Sqdeg_SN${SN}_NoMask_LSSTGpAM_X3sigma_SS${SS}_zLSST_ZBcut${ZBcut}_ThBins9_MRres60arcsec

paramfile=$1
los_start=$2
los_end=$3
./Master_CorrFun_ByParts.sh Sims_Run $paramfile $los_start $los_end


###########################SBATCH --mail-user=bengib@roe.ac.uk 
###########SBATCH -o log_Cosmol${i}.out 
