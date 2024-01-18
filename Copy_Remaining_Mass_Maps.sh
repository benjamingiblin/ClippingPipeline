#!/bin/bash
# Scroll through workers, find mass maps that still exist on them
# and copy them to the head node

MRres="140.64arcs" # 60arcsec
SS=2.816 #1.408 #3.11, 9.33, 18.66
Mask="Mosaic"
Survey="KiDS1000"
SN=(0.27 0.258 0.273 0.254 0.27)
ZBcut=("0.1-0.3" "0.3-0.5" "0.5-0.7" "0.7-0.9" "0.9-1.2")

data_DIR='/data/bengib/Clipping_Pipeline/'
home_DIR='/home/bengib/Clipping_Pipeline/'

for i in `seq 2 80`; do # scroll through workers
    echo "------------------------------- On worker $i -------------------------------"
    printf -v i "%03d" $i
    for cosmol in `seq 0 24`; do
	echo "__________________ cosmol $cosmol ___________________ "
	for j in `seq 0 4`; do
	    jp1=$((j+1))
	    for k in `seq $jp1 4`; do
		if [ "$j" -eq "$k" ]; then ZBlabel=ZBcut${ZBcut[$j]}; else ZBlabel=ZBcut${ZBcut[$j]}_X_ZBcut${ZBcut[$k]}; fi
		echo " xxxxxxxxxxx $ZBlabel xxxxxxxxxxxx "
		
		DIRname1=MRres${MRres}_100Sqdeg_SNCycle_${Mask}_${Survey}GpAM_z${Survey}_${ZBlabel}_Cosmol${cosmol}
		DIRname2=MRres${MRres}_100Sqdeg_NOISE_${Mask}_${Survey}GpAM_z${Survey}_${ZBlabel}
		#ssh worker$i ls $data_DIR/Mass_Recon/$DIRname/*Ekappa.npy
		ssh worker$i rsync -avz $data_DIR/Mass_Recon/$DIRname1/*Ekappa.npy $home_DIR/Mass_Recon/$DIRname1/
		#ssh worker$i rsync -avz $data_DIR/Mass_Recon/$DIRname2/*Ekappa.npy $home_DIR/Mass_Recon/$DIRname2/
		
	    done
	done
    done
done

