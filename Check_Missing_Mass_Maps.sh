#!/bin/bash
# 03/08/2021, B. Giblin
# Check for missing mass maps across all cosmologies, ZBcuts (& combo's), LOS and noise realisations

mock_Type="SLICS" # or "cosmoSLICS"

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


SS=18.66 #3.11, 9.33, 18.66
ZBcut=("0.1-0.3" "0.3-0.5" "0.5-0.7" "0.7-0.9" "0.9-1.2")

# KiDS1000 mocks
Survey="KiDS1000"
SN=(0.27 0.258 0.273 0.254 0.27)

# LSST mocks
#Survey="LSST"
#SN=(0.28 0.28 0.28 0.28 0.28)

# KiDS1000-SLICS and LSST-SLICS have the same missing LOS.
missing_los=(135 140 449 595 596 597 598 601 610 614 735)   

param_dir=/home/bengib/Clipping_SimsLvW/param_files

nn_start=0
nn_end=19  # range of noise realisations to cycle through
nn=0

missing_los=(135 140 449 595 596 597 598 601 610 614 735)
count=0
for i in fid; do
#for i in `seq 0 24`; do   
    echo "--------------------------- COSMOL $i ----------------------------------------"
    for j in `seq 0 4`; do
	jp1=$((j+1));
	echo " xxxxxxxxxxxxxxxxxxxxxxx ${ZBcut[$j]} xxxxxxxxxxxxxxxxxxxxxx";
	for k in `seq $j 4`; do

	    if [ "$j" -eq "$k" ]; then ZBlabel=ZBcut${ZBcut[$j]}; else ZBlabel=ZBcut${ZBcut[$j]}_X_ZBcut${ZBcut[$k]}; fi
	    
	    for los in `seq $los_start $los_end`; do
		#for nn in `seq $nn_start $nn_end`; do
		    # This line skips the missing SLICS los
		    for item in ${missing_los[*]}; do if [[ "$los" -eq "$item" ]]; then los=$((los+1)); fi; done
		

		    if [ "$mock_Type" == "SLICS" ]; then
			# normal map
			#f=Mass_Recon/MRres60arcs_100Sqdeg_SN${SN[$j]}_NoMask_${Survey}GpAM_z${Survey}_${ZBlabel}/SN${SN[$j]}_test.${Survey}GpAM.LOS${los}.SS${SS}.Ekappa.fits

			# a pure noise map
			f=Mass_Recon/MRres60arcs_100Sqdeg_NOISE_NoMask_${Survey}GpAM_z${Survey}_${ZBlabel}/NOISE_test.${Survey}GpAM.LOS${los}.SS${SS}.Ekappa.fits
		    else
			# all the noise realisations
			f=Mass_Recon/MRres60arcs_${IA_Tag}100Sqdeg_SNCycle_NoMask_${Survey}GpAM_z${Survey}_${ZBlabel}_Cosmol${i}/SN${SN[$j]}_test.${Survey}GpAM.LOS${los}n${nn}.SS${SS}.Ekappa.fits

			# one noise real. per LOS
			#f=Mass_Recon/MRres60arcs_${IA_Tag}100Sqdeg_SN${SN[$j]}_NoMask_${Survey}GpAM_z${Survey}_${ZBlabel}_Cosmol${i}/SN${SN[$j]}_test.${Survey}GpAM.LOS${los}.SS${SS}.Ekappa.fits
			
		    fi

		    if [ ! -f $f ]; then

			# NOISE-ONLY
			paramfile1=$param_dir/100Sqdeg_NOISE_NoMask_${Survey}GpAM_X3sigma_SS${SS}_z${Survey}_ZBcut${ZBcut[$j]}_ThBins9_MRres60arcsec
			paramfile2=$param_dir/100Sqdeg_NOISE_NoMask_${Survey}GpAM_X3sigma_SS${SS}_z${Survey}_ZBcut${ZBcut[$k]}_ThBins9_MRres60arcsec                   

			# cosmoSLICS - all the noise real.
			#paramfile1=$param_dir/${IA_Tag}100Sqdeg_CycleSN${SN[$j]}_NoMask_${Survey}GpAM_X3sigma_SS${SS}_z${Survey}_ZBcut${ZBcut[$j]}_ThBins9_Cosmol${i}_MRres60arcsec
			#paramfile2=$param_dir/${IA_Tag}100Sqdeg_CycleSN${SN[$k]}_NoMask_${Survey}GpAM_X3sigma_SS${SS}_z${Survey}_ZBcut${ZBcut[$k]}_ThBins9_Cosmol${i}_MRres60arcsec
			# cosmoSLICS - one noise real. per LOS.
			#paramfile1=$param_dir/100Sqdeg_SN${SN[$j]}_NoMask_${Survey}GpAM_X3sigma_SS${SS}_z${Survey}_ZBcut${ZBcut[$j]}_ThBins9_Cosmol${i}_MRres60arcsec
                        #paramfile2=$param_dir/100Sqdeg_SN${SN[$k]}_NoMask_${Survey}GpAM_X3sigma_SS${SS}_z${Survey}_ZBcut${ZBcut[$k]}_ThBins9_Cosmol${i}_MRres60arcsec
			

			# SLICS
			#paramfile1=$param_dir/100Sqdeg_SN${SN[$j]}_NoMask_${Survey}GpAM_X3sigma_SS${SS}_z${Survey}_ZBcut${ZBcut[$j]}_ThBins9_MRres60arcsec
			#paramfile2=$param_dir/100Sqdeg_SN${SN[$k]}_NoMask_${Survey}GpAM_X3sigma_SS${SS}_z${Survey}_ZBcut${ZBcut[$k]}_ThBins9_MRres60arcsec

			
			# This bit of code re-launches SPECIFIC los and noise realisations
			# (to be used when you're close to getting them all!)
			if [ "$Launch_Specific_LOS" == "True" ]; then
			    echo $i ${SN[$j]} ${ZBcut[$j]} ${SN[$k]} ${ZBcut[$k]} ${SS} ${los}n$nn
			    count=$((count+1))
			    r=$RANDOM
			    #param_files/Change_Specific_Paramfile.sh $paramfile1 $nn $r
			    #param_files/Change_Specific_Paramfile.sh $paramfile2 $nn $r
			    #if [ "$j" -eq "$k" ]; then
				# auto-bin
				#sbatch Launch.sh ${paramfile1}_tmp${r} $los $los
			    #else
				# cross-bin
  				#sbatch Launch_Combine-zbins.sh ${paramfile1}_tmp${r} $los $los ${paramfile2}_tmp${r}
			    #fi
			    
			else
			    echo $i ${SN[$j]} ${ZBcut[$j]} ${SN[$k]} ${ZBcut[$k]} ${SS} ${los} #${los_end}
			    count=$((count+1))
			    if [ "$j" -eq "$k" ]; then
				# auto-bin
				sbatch Launch.sh $paramfile1 $los $los
			    else
				# cross-bin
				sbatch Launch_Combine-zbins.sh $paramfile1 $los $los $paramfile2  
			    fi
			    
			    #break # stop scrolling through LOS now you found a missing one for this cosmol & zbin
			fi
			
		    fi
		#done
	    done
	done
    done
done

echo "$count jobs submitted."
