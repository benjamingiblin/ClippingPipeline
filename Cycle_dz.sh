#!/bin/bash
# compute combo of dzX & dzY results for the redshift cross bins

ZBcut=("0.1-0.3" "0.3-0.5" "0.5-0.7" "0.7-0.9" "0.9-1.2")
SN=(0.27 0.258 0.273 0.254 0.27)
dz=("dz2" "dz3" "dz4" "dz5")

config=config_cross_dz.txt
pdir=/home/bengib/Clipping_Pipeline/param_files/


count=1
for j in `seq 0 4`; do
    jp1=$((j+1))
    for k in `seq $jp1 4`; do
	echo " ------- ZBcut${ZBcut[$j]}_X_ZBcut${ZBcut[$k]} ----------- "
	for sx in `seq 0 3`; do
	    for sy in `seq 0 3`; do  #0 3

		echo "$count ${dz[$sx]}_X_${dz[$sy]}"

		p1=${pdir}/${dz[$sx]}_100Sqdeg_SN${SN[$j]}_Mosaic_KiDS1000GpAM_X3sigma_SS2.816_zKiDS1000_ZBcut${ZBcut[$j]}_ThBins9_Cosmolfid_MRres140.64arcsec
		p2=${pdir}/${dz[$sy]}_100Sqdeg_SN${SN[$k]}_Mosaic_KiDS1000GpAM_X3sigma_SS2.816_zKiDS1000_ZBcut${ZBcut[$k]}_ThBins9_Cosmolfid_MRres140.64arcsec

		# add to config
		#echo $count $p1 74 74 $p2 >> $config
		# or launch!
		sbatch Launch_Combine-zbins.sh $p1 74 74 $p2
		echo "sleeping for 2min..."; sleep 2m
		
		count=$((count+1))
	    done
	done
    done
done
