#!/bin/bash

SN=("0.27" "0.258" "0.273" "0.254" "0.27")
Z=("0.1-0.3" "0.3-0.5" "0.5-0.7" "0.7-0.9" "0.9-1.2")

for i in `seq 0 4`; do
    iplus1=$((i+1))
    for j in `seq $iplus1 4`; do

	#sbatch Launch.sh /home/bengib/Clipping_Pipeline/param_files/IA1.0_100Sqdeg_SN${SN[$i]}_NoMask_KiDS1000GpAM_X3sigma_SS2.816_zKiDS1000_ZBcut${Z[$i]}_ThBins9_Cosmolfid_MRres140.64arcsec 1 50 /home/bengib/Clipping_Pipeline/param_files/IA1.0_100Sqdeg_SN${SN[$j]}_NoMask_KiDS1000GpAM_X3sigma_SS2.816_zKiDS1000_ZBcut${Z[$j]}_ThBins9_Cosmolfid_MRres140.64arcsec 

	echo "${SN[$i]} ${Z[$i]} ----- ${SN[$j]} ${Z[$j]}"
    done
done
