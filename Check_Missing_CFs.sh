#!/bin/bash
# 03/08/2021, B. Giblin
# Check for missing CFs cross all cosmologies, ZBcuts (& combo's), LOS and noise realisations

mock_Type="cosmoSLICS" # or "cosmoSLICS"

Launch_Specific_LOS="True" #"True" / "False"

# If False, it runs everything from the missing LOS to the end (incl all noise realisations
# If True, it runs launches jobs with specific LOS and noise realisations
# It's only worth doing the latter when you've got ALMOST all CFs
# Otherwise it'll launch a ridiculous number of (admittedly fast) jobs.


if [ "$mock_Type" == "SLICS" ]; then
    los_start=74
    los_end=799 #699 #799 for KiDS1000, 699 for LSST
else
    los_start=1
    los_end=50
fi


Check_Sys="BaryOFF" # "IA", "BaryON", "BaryOFF", "dz" or "None"

if [ "$Check_Sys" == "IA" ]; then
    echo " !!! SEARCHING FOR THE IA MEASUREMENTS !!! "
    ia=$2
    Sys_Tag="IA${ia}_"
    Sys_Tag1=$Sys_Tag
    Sys_Tag2=$Sys_Tag
    los_start=1
    los_end=50

elif [[ "$Check_Sys" == *"Bary"* ]]; then
    echo " !!! SEARCHING FOR THE ${Check_Sys} MEASUREMENTS !!! "
    Sys_Tag="${Check_Sys}_"
    Sys_Tag1=$Sys_Tag
    Sys_Tag2=$Sys_Tag
    los_start=0
    los_end=9

elif [ "$Check_Sys" == "dz" ]; then
    echo " !!! SEARCHING FOR THE ${Check_Sys} MEASUREMENTS !!! "
    dz_vals=("dz2" "dz3" "dz4" "dz5") # option to cycle through combos (see loops below)
    dz=$2 # or use a single dz value read in (MUST COMMENT OUT LOOPS BELOW)
    Sys_Tag="dz${dz}_"
    los_start=74
    los_end=74
fi


MRres="140.64arcs" # 60arcsec
SS=2.816 #3.11, 9.33, 18.66
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
TCdir=/home/bengib/Clipping_Pipeline/Tree_Correlation_Function/

nn_start=0
nn_end=9  # range of noise realisations to cycle through
nn=0

RR_start=1
RR_end=18
RR=2       # region; used if working with mosaic mocks

index=$1
config1=config_auto_${index}.txt
config2=config_cross_${index}.txt
count=0
for file_type in "SS${SS}.rCLIP_X3sigma.CorrFun.asc" "ORIG.CorrFun.asc"; do
    echo " FILE_TYPE: $file_type "
for i in fid; do
#for i in `seq 0 24`; do   
    echo "--------------------------- COSMOL $i ----------------------------------------"
    for j in `seq 0 4`; do
	jp1=$((j+1));
	echo " xxxxxxxxxxxxxxxxxxxxxxx ${ZBcut[$j]} xxxxxxxxxxxxxxxxxxxxxx";
	for k in `seq $jp1 4`; do  # $j 4

	    if [ "$j" -eq "$k" ]; then ZBlabel=ZBcut${ZBcut[$j]}; else ZBlabel=ZBcut${ZBcut[$j]}_X_ZBcut${ZBcut[$k]}; fi
	    
	    for los in `seq $los_start $los_end`; do
		#for nn in `seq $nn_start $nn_end`; doß
		for RR in `seq $RR_start $RR_end`; do 


                    # CYCLE THROUGH dz COMBOS...:                                                                                    
                    #for dz1 in `seq 0 3`; do for dz2 in `seq 0 3`; do
                    #        if [ "$dz1" -eq "$dz2" ]; then
                    #            Sys_Tag=${dz_vals[$dz1]}_; Sys_Tag1=${dz_vals[$dz1]}_; Sys_Tag2=$Sys_Tag1;
                    #        else
                    #            Sys_Tag=${dz_vals[$dz1]}_X_${dz_vals[$dz2]}_;
                    #            Sys_Tag1=${dz_vals[$dz1]}_; Sys_Tag2=${dz_vals[$dz2]}_;
                    #        fi

		    
		    # This line skips the missing SLICS los
		    for item in ${missing_los[*]}; do if [[ "$los" -eq "$item" ]]; then los=$((los+1)); fi; done
		

		    if [ "$mock_Type" == "SLICS" ]; then
			# normal CF
			#f=${TCdir}/MRres${MRres}_${Sys_Tag}100Sqdeg_SN${SN[$j]}_${Mask}_${Survey}GpAM_z${Survey}_${ZBlabel}/ThBins9/SN${SN[$j]}_${filetag}.${Survey}GpAM.LOS${los}R${RR}.${file_type}

			# noise-real. SLICS map:
			f=${TCdir}/MRres${MRres}_100Sqdeg_SNCycle_${Mask}_${Survey}GpAM_z${Survey}_${ZBlabel}/ThBins9/SN${SN[$j]}_${filetag}.${Survey}GpAM.LOS${los}R${RR}n${nn}.${file_type}

			# a pure noise CF
			#f=${TCdir}/MRres${MRres}_${Sys_Tag}100Sqdeg_NOISE_${Mask}_${Survey}GpAM_z${Survey}_${ZBlabel}/ThBins9/NOISE_${filetag}.${Survey}GpAM.LOS${los}R${RR}.${file_type}
		    else
			# all the noise realisations
			#f=${TCdir}/MRres${MRres}_${Sys_Tag}100Sqdeg_SN${SN[$j]}_${Mask}_${Survey}GpAM_z${Survey}_${ZBlabel}_Cosmol${i}/ThBins9/SN${SN[$j]}_${filetag}.${Survey}GpAM.LOS${los}n${nn}.${file_type}

			# one noise real. per LOS
			f=${TCdir}/MRres${MRres}_${Sys_Tag}100Sqdeg_SN${SN[$j]}_${Mask}_${Survey}GpAM_z${Survey}_${ZBlabel}_Cosmol${i}/ThBins9/SN${SN[$j]}_${filetag}.${Survey}GpAM.LOS${los}R${RR}.${file_type}
			
		    fi

		    if [ ! -f $f ]; then

			# NOISE-ONLY
			#paramfile1=$param_dir/100Sqdeg_NOISE_${Mask}_${Survey}GpAM_X3sigma_SS${SS}_z${Survey}_ZBcut${ZBcut[$j]}_ThBins9_MRres${MRres}ec
			#paramfile2=$param_dir/100Sqdeg_NOISE_${Mask}_${Survey}GpAM_X3sigma_SS${SS}_z${Survey}_ZBcut${ZBcut[$k]}_ThBins9_MRres${MRres}ec                   

			# cosmoSLICS - all the noise real.
			#paramfile1=$param_dir/${Sys_Tag1}100Sqdeg_CycleSN${SN[$j]}_${Mask}_${Survey}GpAM_X3sigma_SS${SS}_z${Survey}_ZBcut${ZBcut[$j]}_ThBins9_Cosmol${i}_MRres${MRres}ec
			#paramfile2=$param_dir/${Sys_Tag2}100Sqdeg_CycleSN${SN[$k]}_${Mask}_${Survey}GpAM_X3sigma_SS${SS}_z${Survey}_ZBcut${ZBcut[$k]}_ThBins9_Cosmol${i}_MRres${MRres}ec
			# cosmoSLICS - one noise real. per LOS.
			paramfile1=$param_dir/${Sys_Tag1}100Sqdeg_SN${SN[$j]}_${Mask}_${Survey}GpAM_X3sigma_SS${SS}_z${Survey}_ZBcut${ZBcut[$j]}_ThBins9_Cosmol${i}_MRres${MRres}ec
                        paramfile2=$param_dir/${Sys_Tag2}100Sqdeg_SN${SN[$k]}_${Mask}_${Survey}GpAM_X3sigma_SS${SS}_z${Survey}_ZBcut${ZBcut[$k]}_ThBins9_Cosmol${i}_MRres${MRres}ec
			

			# SLICS - one noise real. per LOS. 
			#paramfile1=$param_dir/100Sqdeg_SN${SN[$j]}_${Mask}_${Survey}GpAM_X3sigma_SS${SS}_z${Survey}_ZBcut${ZBcut[$j]}_ThBins9_MRres${MRres}ec
			#paramfile2=$param_dir/100Sqdeg_SN${SN[$k]}_${Mask}_${Survey}GpAM_X3sigma_SS${SS}_z${Survey}_ZBcut${ZBcut[$k]}_ThBins9_MRres${MRres}ec
			# SLICS - all the noise real.
			#paramfile1=$param_dir/100Sqdeg_CycleSN${SN[$j]}_${Mask}_${Survey}GpAM_X3sigma_SS${SS}_z${Survey}_ZBcut${ZBcut[$j]}_ThBins9_MRres${MRres}ec
			#paramfile2=$param_dir/100Sqdeg_CycleSN${SN[$k]}_${Mask}_${Survey}GpAM_X3sigma_SS${SS}_z${Survey}_ZBcut${ZBcut[$k]}_ThBins9_MRres${MRres}ec
			#ls $paramfile1

			
			# This bit of code re-launches SPECIFIC los and noise realisations
			# (to be used when you're close to getting them all!)
			if [ "$Launch_Specific_LOS" == "True" ]; then
			    echo $f
			    echo $i ${SN[$j]} ${ZBcut[$j]} ${SN[$k]} ${ZBcut[$k]} ${SS} ${los}R${RR} #n$nn
			    count=$((count+1))
			    r=$RANDOM			    
			    param_files/Change_Specific_Paramfile.sh $paramfile1 $RR $r $nn
			    param_files/Change_Specific_Paramfile.sh $paramfile2 $RR $r $nn
			    if [ "${ZBcut[$j]}" == "0.1-1.2" ]; then
				# we need to alter paramfiles for the 5 indiv. zbins:
				ZB_array=("0.1-0.3" "0.3-0.5" "0.5-0.7" "0.7-0.9" "0.9-1.2")
				SN_array=(0.27 0.258 0.273 0.254 0.27)
				for z_idx in `seq 0 4`; do
				    tmp_paramfile="${paramfile1/ZBcut0.1-1.2/ZBcut${ZB_array[$z_idx]}}"   # replace the ZBcut label
				    tmp_paramfile="${tmp_paramfile/SN${SN[$j]}_/SN${SN_array[$z_idx]}_}"     # replace the SN label 
				    param_files/Change_Specific_Paramfile.sh $tmp_paramfile $RR $r
				done
			    fi
				    
			    if [ "$j" -eq "$k" ]; then
				# auto-bin
				#sbatch Launch.sh ${paramfile1}_tmp${r} $los $los
				echo "$count ${paramfile1}_tmp${r} $los $los" >> $config1
			    else
				# cross-bin
				#sbatch Launch.sh ${paramfile1}_tmp${r} $los $los # if the *.asc cat doesnt exist
				#sbatch Launch.sh ${paramfile2}_tmp${r} $los $los # if the *.asc cat doesnt exist
				#echo "$count ${paramfile1}_tmp${r} ${paramfile2}_tmp${r} $los $los $RR $RR $nn $nn" >> $config2
				# or this one, to calc the CF itself:
				sbatch Launch_Calc_CrossCorr_Cycle_z_And_LOS.sh ${paramfile1}_tmp${r} ${paramfile2}_tmp${r} $los $los $RR $RR
				sleep 2s
			    fi
			    
			else
			    echo $f
			    echo $i ${SN[$j]} ${ZBcut[$j]} ${SN[$k]} ${ZBcut[$k]} ${SS} ${los}R${RR} #${los_end}
			    count=$((count+1))
			    if [ "$j" -eq "$k" ]; then
				# auto-bin
				#sbatch Launch.sh $paramfile1 $los $los
				#sbatch Launch_Calc_CrossCorr_Cycle_z_And_LOS.sh $paramfile1 $paramfile2 $los $los
				echo "$count $paramfile1 $paramfile2 $los $los $R $R" >> $config1
			    else
				# cross-bin
			    	#sbatch Launch.sh $paramfile1 $los $los
				#sbatch Launch.sh $paramfile2 $los $los
				sbatch Launch_Calc_CrossCorr_Cycle_z_And_LOS.sh $paramfile1 $paramfile2 $los_start $los_end 1 18
				echo "$count $paramfile1 $paramfile2 $los_start $los_end 1 18" >> $config2
				#count=$((count+1))
				#echo "$count $paramfile2 $los $los" >> $config2
			    fi
			    
			#    break # stop scrolling through LOS now you found a missing one for this cosmol & zbin
			fi
			
		    fi
		    #done; done # CLOSE dz COMBO LOOPS
		done # RR loop
		#done  # nn loop
		    #python Correlation_Function/plot_CorrFun.py Sims_Run $paramfile1 $los_start $los_end $paramfile2
		    #break
	    done
	done
    done
done
done

echo "$count jobs submitted."
