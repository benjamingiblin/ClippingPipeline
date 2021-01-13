#!/bin/bash

overall_DIR=/home/bengib/Clipping_SimsLvW/

start=1
end=157
Noise_start=0 # Noise real you started on
Noise_end=60  # ... you finished on
Noise="$((Noise_end - Noise_start +1))" # Number Noise real in this range PER CATALOGUE
NLOS="$((5 * Noise))"        # Number of noise realisations that should exist IN TOTAL
DIRname='36Sqdeg_SNCycle_NoMask_3.32GpAM_zKiDS_ZBcut0.5-0.9'    # DIRname structure excluding cosmol
DATE=$(date +\%Y-\%m-\%d)
SLURM_DIR=SLURM_LOGS/SLURM_LOGS_Noise0-60_NewKappac      

# Find Missing Files or Re-run the whole thing?
FindFiles="Y"
ReRun="N"

for i in `seq $start $end`; do
    if [ ! -d "$overall_DIR/Tree_Correlation_Function/${DIRname}_Cosmol${i}/ThBins9/NLOS$NLOS" ]; then
	string=$(grep -w "Subdir: ${DIRname}_Cosmol${i}" $overall_DIR/$SLURM_DIR/slurm*)
	string=$(echo ${string%: Output* })
	slurmfile=`basename $string`                     # extract name of slurm file
	worker=$(grep "worker" $overall_DIR/$SLURM_DIR/$slurmfile)  # extract name of worker
	echo "Cosmology $i, slurm-file $slurmfile, did not complete. Removing files from $worker"
        echo "Cosmology $i, slurm-file $slurmfile, did not complete. Removing files from $worker" >> $overall_DIR/Tree_Correlation_Function/Find_Failures_Log$DATE.txt

	# Find which output CF files do not exist...?
	if [ "$FindFiles" == "Y" ]; then
	    echo "Looking for missing files...."
	    for l in `seq 0 4`; do
		if [[ "$DIRname" == *"SNCycle"* ]]; then           # If you are searching for noise realisations, output filename is different
#		    nminus1="$((Noise - 1))"
		    for n in `seq $Noise_start $Noise_end`; do
			if [ ! -f "$overall_DIR/Tree_Correlation_Function/${DIRname}_Cosmol${i}/ThBins9/SN0.28_test.3.32GpAM.LOS${l}n${n}.SS112.rCLIP_X3sigma.CorrFun.asc" ]; then
			    echo "        Cosmology ${i}: ${l}n${n} is missing."
			    echo "        Cosmology ${i}: ${l}n${n} is missing." >> $overall_DIR/Tree_Correlation_Function/Find_Failures_Log$DATE.txt
			fi
		    done
		else
		    if [ ! -f "$overall_DIR/Tree_Correlation_Function/${DIRname}_Cosmol${i}/ThBins9/SN0.28_test.3.32GpAM.LOS${l}.SS112.rCLIP_X3sigma.CorrFun.asc" ]; then
			echo "        Cosmology ${i}: ${l} is missing."
			echo "        Cosmology ${i}: ${l} is missing." >> $overall_DIR/Tree_Correlation_Function/Find_Failures_Log$DATE.txt
		    fi
		fi
	    done
	fi
	
	
	# Rerun the whole thing...?
	if [ "$ReRun" == "Y" ]; then sbatch $overall_DIR/Launch.sh $i; fi

	
	# Conservative cleanse of worker - RUN THIS IF YOU HAVE JOBS RUNNING
	#ssh $worker rm -rf /data/bengib/Clipping_SimsLvW/*/${DIRname}_Cosmol${i}
	# !!! AGGRESSIVE !!! cleanse of worker - DON'T DO THIS IF YOU HAVE JOBS RUNNING
	#ssh $worker rm -rf /data/bengib/Clipping_SimsLvW/ 

    fi
done
TIME=$(date +\%H:\%M)
echo "----------------------------- RAN AT $TIME Cosmol: $start to $end------------------------------" >> $overall_DIR/Tree_Correlation_Function/Find_Failures_Log$DATE.txt





