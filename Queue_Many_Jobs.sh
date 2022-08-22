#!/bin/bash

# 12/11/2020, B. Giblin
# Launch many clipping/mass mapping jobs, possibly varying cosmology, ZBcut and hence sigma_e.
# Important: Launch.sh must be set up to receive SN & ZBcut params.
# Also make sure the params file in Launch.sh is set to the desired one.

ZBcut=("0.1-0.3" "0.3-0.5" "0.5-0.7" "0.7-0.9" "0.9-1.2")
SN=(0.27 0.258 0.273 0.254 0.27)
SS=18.66

echo "Sleeping for 2hrs"
#sleep 2h

echo "Removing old files"
#Tree_Correlation_Function/Clean_Workers.sh
#sleep 15m


# launch jobs for each cosmology and redshift cut:
for i in `seq 0 24`; do
    for j in `seq 0 4`; do
	jp1=$((j+1))
	#for k in `seq $jp1 4`; do
	for k in `seq $j $j`; do   
	    p1=/home/bengib/Clipping_SimsLvW/param_files/100Sqdeg_CycleSN${SN[$j]}_NoMask_KiDS1000GpAM_X3sigma_SS${SS}_zKiDS1000_ZBcut${ZBcut[$j]}_ThBins9_Cosmol${i}_MRres60arcsec
	    p2=/home/bengib/Clipping_SimsLvW/param_files/100Sqdeg_CycleSN${SN[$k]}_NoMask_KiDS1000GpAM_X3sigma_SS${SS}_zKiDS1000_ZBcut${ZBcut[$k]}_ThBins9_Cosmol${i}_MRres60arcsec

	    #ls $p1 $p2
	    echo "Running cosmol $i, ZBcut ${ZBcut[$j]}_X_${ZBcut[$k]}, sigma_e ${SN[$j]}"
	    #sbatch Launch_Combine-zbins.sh $p1 1 50 $p2
	    sbatch Launch.sh $p1 1 50 

	    if [ "$i" == "0" ]; then
		echo "RUNNING FID AND NOISE"
		# run fid and noise
		pfid1=/home/bengib/Clipping_SimsLvW/param_files/100Sqdeg_CycleSN${SN[$j]}_NoMask_KiDS1000GpAM_X3sigma_SS${SS}_zKiDS1000_ZBcut${ZBcut[$j]}_ThBins9_Cosmolfid_MRres60arcsec
		pfid2=/home/bengib/Clipping_SimsLvW/param_files/100Sqdeg_CycleSN${SN[$k]}_NoMask_KiDS1000GpAM_X3sigma_SS${SS}_zKiDS1000_ZBcut${ZBcut[$k]}_ThBins9_Cosmolfid_MRres60arcsec
		#ls $pfid1 $pfid2
		#sbatch Launch_Combine-zbins.sh $pfid1 1 50 $pfid2
		sbatch Launch.sh $pfid1 1 50
		
		pn1=/home/bengib/Clipping_SimsLvW/param_files/100Sqdeg_NOISE_NoMask_KiDS1000GpAM_X3sigma_SS${SS}_zKiDS1000_ZBcut${ZBcut[$j]}_ThBins9_MRres60arcsec
		pn2=/home/bengib/Clipping_SimsLvW/param_files/100Sqdeg_NOISE_NoMask_KiDS1000GpAM_X3sigma_SS${SS}_zKiDS1000_ZBcut${ZBcut[$k]}_ThBins9_MRres60arcsec
		#ls $pn1 $pn2
		#sbatch Launch_Combine-zbins.sh $pn1 74 799 $pn2
		#sbatch Launch.sh $pn1 74 799
	    fi
	    
	    
	done
    done
done


