#!/bin/bash
# 03/08/2021, B. Giblin
# Check for missing mass maps across all cosmologies, ZBcuts (& combo's), LOS and noise realisations

mock_Type="cosmoSLICS" # or "cosmoSLICS"

Launch_Specific_LOS="False" #"True" / "False"

# If False, it runs everything from the missing LOS to the end (incl all noise realisations
# If True, it runs launches jobs with specific LOS and noise realisations
# It's only worth doing the latter when you've got ALMOST all maps
# Otherwise it'll launch a ridiculous number of (admittedly fast) jobs.


if [ "$mock_Type" == "SLICS" ]; then
    los_start=74
    los_end=799 #699 #799 for KiDS1000, 699 for LSST
else
    los_start=1
    los_end=50
fi


Check_IA="False" # Look for the maps contaminated by IAs
if [ "$Check_IA" == "True" ]; then
    echo " !!! SEARCHING FOR THE IA-CONTAMINATED MEASUREMENTS !!! "
    IA_Tag="IA0.0_" # (else it's blank)
fi


MRres="140.64arcs" # 60arcsec
SS=2.816 #1.408, 2.816, 5.631 #3.11, 9.33, 18.66
Mask="Mosaic"
if [ "$Mask" == "Mosaic" ]; then
    filetag="Mosaic"
    missing_los=(198 199) # Mosaic mocks have different missing LOS.
    if [ "$mock_Type" == "SLICS" ]; then los_end=292; fi
else
    filetag="test"
    missing_los=(135 140 449 595 596 597 598 601 610 614 735)
fi

ZBcut=("0.1-0.3" "0.3-0.5" "0.5-0.7" "0.7-0.9" "0.9-1.2")
#ZBcut=("0.1-1.2")

# KiDS1000 mocks
Survey="KiDS1000"
SN=(0.27 0.258 0.273 0.254 0.27)
#SN=(0.265)

# LSST mocks
#Survey="LSST"
#SN=(0.28 0.28 0.28 0.28 0.28)

param_dir=/home/bengib/Clipping_Pipeline/param_files

nn_start=0
nn_end=9  # range of noise realisations to cycle through
nn=0

RR_start=1
RR_end=18
RR=1       # region; used if working with mosaic mocks

index=$1
# make a files to submit job arrays
config1=config_auto_${index}.txt
config2=config_cross_${index}.txt
echo "ArrayTaskID PARAMFILE LOS_START LOS_END" > $config1
echo "ArrayTaskID PARAMFILE1 LOS_START LOS_END PARAMFILE2" > $config2

count=0
#for i in fid; do
for i in `seq $index $index`; do   
    echo "--------------------------- COSMOL $i ----------------------------------------"
    for j in `seq 0 4`; do
	jp1=$((j+1));
	echo " xxxxxxxxxxxxxxxxxxxxxxx ${ZBcut[$j]} xxxxxxxxxxxxxxxxxxxxxx";
	for k in `seq $jp1 4`; do  # $j 4

	    if [ "$j" -eq "$k" ]; then ZBlabel=ZBcut${ZBcut[$j]}; else ZBlabel=ZBcut${ZBcut[$j]}_X_ZBcut${ZBcut[$k]}; fi
	    
	    for los in `seq $los_start $los_end`; do
		#for RR in `seq $RR_start $RR_end`; do
		#for nn in `seq $nn_start $nn_end`; do
		
		    # This line skips the missing SLICS los
		    for item in ${missing_los[*]}; do if [[ "$los" -eq "$item" ]]; then los=$((los+1)); fi; done
		

		    if [ "$mock_Type" == "SLICS" ]; then
			# normal map
			f=Mass_Recon/MRres${MRres}_100Sqdeg_SN${SN[$j]}_${Mask}_${Survey}GpAM_z${Survey}_${ZBlabel}/SN${SN[$j]}_${filetag}.${Survey}GpAM.LOS${los}R${RR}.SS${SS}.Ekappa.npy

			# a pure noise map
			#f=Mass_Recon/MRres${MRres}_100Sqdeg_NOISE_${Mask}_${Survey}GpAM_z${Survey}_${ZBlabel}/NOISE_${filetag}.${Survey}GpAM.LOS${los}R${RR}.SS${SS}.Ekappa.npy
		    else
			# all the noise realisations
			f=Mass_Recon/MRres${MRres}_${IA_Tag}100Sqdeg_SNCycle_${Mask}_${Survey}GpAM_z${Survey}_${ZBlabel}_Cosmol${i}/SN${SN[$j]}_${filetag}.${Survey}GpAM.LOS${los}R${RR}n${nn}.SS${SS}.Ekappa.npy

			# one noise real. per LOS
			#f=Mass_Recon/MRres${MRres}_${IA_Tag}100Sqdeg_SN${SN[$j]}_${Mask}_${Survey}GpAM_z${Survey}_${ZBlabel}_Cosmol${i}/SN${SN[$j]}_${filetag}.${Survey}GpAM.LOS${los}R${RR}.SS${SS}.Ekappa.npy
			
		    fi

		    if [ ! -f $f ]; then

			# NOISE-ONLY
			#paramfile1=$param_dir/100Sqdeg_NOISE_${Mask}_${Survey}GpAM_X3sigma_SS${SS}_z${Survey}_ZBcut${ZBcut[$j]}_ThBins9_MRres${MRres}ec
			#paramfile2=$param_dir/100Sqdeg_NOISE_${Mask}_${Survey}GpAM_X3sigma_SS${SS}_z${Survey}_ZBcut${ZBcut[$k]}_ThBins9_MRres${MRres}ec                   

			# cosmoSLICS - all the noise real.
			paramfile1=$param_dir/${IA_Tag}100Sqdeg_CycleSN${SN[$j]}_${Mask}_${Survey}GpAM_X3sigma_SS${SS}_z${Survey}_ZBcut${ZBcut[$j]}_ThBins9_Cosmol${i}_MRres${MRres}ec
			paramfile2=$param_dir/${IA_Tag}100Sqdeg_CycleSN${SN[$k]}_${Mask}_${Survey}GpAM_X3sigma_SS${SS}_z${Survey}_ZBcut${ZBcut[$k]}_ThBins9_Cosmol${i}_MRres${MRres}ec
			# cosmoSLICS - one noise real. per LOS.
			#paramfile1=$param_dir/100Sqdeg_SN${SN[$j]}_${Mask}_${Survey}GpAM_X3sigma_SS${SS}_z${Survey}_ZBcut${ZBcut[$j]}_ThBins9_Cosmol${i}_MRres${MRres}ec
                        #paramfile2=$param_dir/100Sqdeg_SN${SN[$k]}_${Mask}_${Survey}GpAM_X3sigma_SS${SS}_z${Survey}_ZBcut${ZBcut[$k]}_ThBins9_Cosmol${i}_MRres${MRres}ec
			

			# SLICS
			#paramfile1=$param_dir/100Sqdeg_SN${SN[$j]}_${Mask}_${Survey}GpAM_X3sigma_SS${SS}_z${Survey}_ZBcut${ZBcut[$j]}_ThBins9_MRres${MRres}ec
			#paramfile2=$param_dir/100Sqdeg_SN${SN[$k]}_${Mask}_${Survey}GpAM_X3sigma_SS${SS}_z${Survey}_ZBcut${ZBcut[$k]}_ThBins9_MRres${MRres}ec
			#ls $paramfile1

			
			# This bit of code re-launches SPECIFIC los and noise realisations
			# (to be used when you're close to getting them all!)
			if [ "$Launch_Specific_LOS" == "True" ]; then
			    echo $i ${SN[$j]} ${ZBcut[$j]} ${SN[$k]} ${ZBcut[$k]} ${SS} ${los}R${RR}n$nn
			    count=$((count+1))
			    r=$RANDOM
			    #param_files/Change_Specific_Paramfile.sh $paramfile1 $RR $r
			    #param_files/Change_Specific_Paramfile.sh $paramfile2 $RR $r
			    if [ "$j" -eq "$k" ]; then
				# auto-bin
				#sbatch Launch.sh ${paramfile1}_tmp${r} $los $los
				echo "$count ${paramfile1}_tmp${r} $los $los" >> $config1
			    else
				# cross-bin
  				#sbatch Launch_Combine-zbins.sh ${paramfile1}_tmp${r} $los $los ${paramfile2}_tmp${r}
				echo "$count ${paramfile1}_tmp${r} $los $los ${paramfile2}_tmp${r}" >> $config2
			    fi
			    #break
			    
			else
			    echo $i ${SN[$j]} ${ZBcut[$j]} ${SN[$k]} ${ZBcut[$k]} ${SS} ${los} #R${RR}n${nn} #_start} ${los_end}
			    count=$((count+1))
			    if [ "$j" -eq "$k" ]; then
				# auto-bin
				echo "$count $paramfile1 $los_start $los_end" >> $config1
				#sbatch Launch.sh $paramfile1 $los_start $los_end
			    else
				# cross-bin
				echo "$count $paramfile1 $los_start $los_end $paramfile2" >> $config2
				#sbatch Launch_Combine-zbins.sh $paramfile1 $los_start $los_end $paramfile2  
			    fi

			    #echo "waiting 20 min to give jobs a head start"
			    #sleep 20m
			    #break # stop scrolling through LOS now you found a missing one for this cosmol & zbin
			fi
			
		    fi
		#done
		#done
	    done
	done
    done
done

echo "$count jobs submitted."

