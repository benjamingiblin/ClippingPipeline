#!/bin/bash
# 12/10/2023: Make a list of inputs (paramfiles, los) and save to a file
# this can then be run by Launch_JobArray.sh

DIR=/home/bengib/Clipping_Pipeline/param_files
SN=(0.27 0.258 0.273 0.254 0.27)
ZBcut=("0.1-0.3" "0.3-0.5" "0.5-0.7" "0.7-0.9" "0.9-1.2")

SS=2.816 #1.408 #5.631
los_start=1
los_end=50


config_file=config3.txt
echo "ArrayTaskID PARAMFILE LOS_START LOS_END" > $config_file     # single paramfile
#echo "ArrayTaskID PARAMFILE1 LOS_START LOS_END PARAMFILE2" > $config_file  # double paramfile

count=1
#for i in fid; do
for i in `seq 15 24`; do
    for j in `seq 0 4`; do
	jp1=$((j+1))
	p1=$DIR/100Sqdeg_CycleSN${SN[$j]}_Mosaic_KiDS1000GpAM_X3sigma_SS${SS}_zKiDS1000_ZBcut${ZBcut[$j]}_ThBins9_Cosmol${i}_MRres140.64arcsec
	for k in `seq $j $j`; do
	    p2=$DIR/100Sqdeg_CycleSN${SN[$k]}_Mosaic_KiDS1000GpAM_X3sigma_SS${SS}_zKiDS1000_ZBcut${ZBcut[$k]}_ThBins9_Cosmol${i}_MRres140.64arcsec
	    # For cosmoSLICS:
	    for l in `seq 0 4`; do
		los_start=$((1+l*10))
		los_end=$((los_start+9))
		echo "$count $p1 $los_start $los_end" >> $config_file
		count=$((count+1))
	    done
	    
	    #echo "$count $p1 $los_start $los_end" >> $config_file
	    #echo "$count $p1 $los_start $los_end $p2" >> $config_file
	    #count=$((count+1))
	    
	    # For SLICS:
	    #for l in `seq 0 43`; do
		#los_start=$((74+l*5))
		#los_end=$((los_start+4))
		#echo "$count $p1 $los_start $los_end" >> $config_file
		#echo "$count $p1 $los_start $los_end $p2" >> $config_file
		#count=$((count+1))
	    #done
	    
	done
    done
done

