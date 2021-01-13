#!/bin/bash

# 27/06/2019 - change multiple param files in a loop

ZBcut=("0.1-0.3" "0.3-0.5" "0.5-0.7" "0.7-0.9" "0.9-1.2")
SN=(0.27 0.258 0.273 0.254 0.27)

SS=1
for i in `seq 0 24`; do       # scroll through cosmologies
    for j in in `seq 0 4`; do
	
	cp 100Sqdeg_SN${SN[$j]}_NoMask_KiDS1000GpAM_X3sigma_SS${SS}_zKiDS1000_ZBcut${ZBcut[$j]}_ThBins9_Cosmolfid_MRres60arcsec 100Sqdeg_SN${SN[$j]}_NoMask_KiDS1000GpAM_X3sigma_SS${SS}_zKiDS1000_ZBcut${ZBcut[$j]}_ThBins9_Cosmol${i}_MRres60arcsec

   
	find 100Sqdeg_SN${SN[$j]}_NoMask_KiDS1000GpAM_X3sigma_SS${SS}_zKiDS1000_ZBcut${ZBcut[$j]}_ThBins9_Cosmol${i}_MRres60arcsec -exec sed -i "s#fid *#$i #" {} +
	
    done
done


