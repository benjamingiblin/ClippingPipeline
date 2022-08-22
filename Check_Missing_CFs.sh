#!/bin/bash
# 03/08/2021, B. Giblin
# Check for missing mass maps across all cosmologies, ZBcuts (& combo's), LOS and noise realisations

mock_Type="cosmoSLICS" # or "cosmoSLICS"

if [ "$mock_Type" == "SLICS" ]; then
    los_start=74
    los_end=699 #799 for KiDS1000, 699 for LSST
else
    los_start=1
    los_end=50
fi


Check_IA="True" # Look for the maps contaminated by IAs
if [ "$Check_IA" == "True" ]; then
    echo " !!! SEARCHING FOR THE IA-CONTAMINATED MEASUREMENTS !!! "
    IA_Tag="IA0.0_" # (else it's blank)
fi



SS=9.33
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
		
		    inDIR=Tree_Correlation_Function/MRres60arcs_${IA_Tag}100Sqdeg_SN${SN[$j]}_NoMask_${Survey}GpAM_z${Survey}_${ZBlabel}
		    if [ "$mock_Type" == "cosmoSLICS" ]; then
			cosmol_label="_Cosmol${i}"
		    else
			cosmol_label=""
		    fi
		    
		    inDIR=${inDIR}${cosmol_label}
                    paramfile1=$param_dir/${IA_Tag}100Sqdeg_SN${SN[$j]}_NoMask_${Survey}GpAM_X3sigma_SS${SS}_z${Survey}_ZBcut$\
{ZBcut[$j]}_ThBins9${cosmol_label}_MRres60arcsec
                    paramfile2=$param_dir/${IA_Tag}100Sqdeg_SN${SN[$k]}_NoMask_${Survey}GpAM_X3sigma_SS${SS}_z${Survey}_ZBcut$\
{ZBcut[$k]}_ThBins9${cosmol_label}_MRres60arcsec
		    
		    # Unclipped
		    f=${inDIR}/ThBins9/SN${SN[$j]}_test.${Survey}GpAM.LOS${los}.ORIG.CorrFun.asc
		    # Clipped
		    #f=${inDIR}/ThBins9/SN${SN[$j]}_test.${Survey}GpAM.LOS${los}.SS${SS}.rCLIP_X3sigma.CorrFun.asc

		    
		    if [ ! -f $f ]; then
			echo $i ${SN[$j]} ${ZBcut[$j]} ${SN[$k]} ${ZBcut[$k]} ${SS} $los $los
			#echo $f
			count=$((count+1))
			#sbatch Launch_Calc_CrossCorr_Cycle_z_And_LOS.sh $paramfile1 $paramfile2 $los $los
			break
		    fi

		    # Use these if you're avg'ing all CFs once you've checked all LOS exist.
		    #sbatch Launch_Calc_CrossCorr_Cycle_z_And_LOS.sh $paramfile1 $paramfile2 $los_start $los_end
		    #python Correlation_Function/plot_CorrFun.py Sims_Run $paramfile1 $los_start $los_end $paramfile2
		    #break
		    
		
		#done
	    done
	    
	done
    done
done

echo "$count jobs submitted."
