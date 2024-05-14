#!/bin/bash
# 05/11/15, B. M. Giblin, PhD student Edinburgh

# The Master Script loops over no. of LOS/Fields running 
# SkelePipeline.sh for each one <-- that contains the breakdown of 
# the clipping pipeline. 



pipeline_DIR='/home/bengib/Clipping_Pipeline/'
data_DIR='/data/bengib/Clipping_Pipeline/'

if [[ "$5" == *"param_files"* ]]; then
    # 2 input paramfiles inputted, this means we're combing shear cats for different zbins
    # assemble the new (combined-redshift) DIRname & filter inputs:
    source $pipeline_DIR/ShowSumClass/Assemble_Combine-zbin-DIRname.sh $1 $2 $3 $4 $5
    paramfile=$paramfile1
else
    # 1 paramfile inputted, just working with shear cat for single zbin
    source $pipeline_DIR/ShowSumClass/FilterInputArgs.sh $1 $2 $3 $4 $5
fi


# The locking mechanism to prevent multiple mass mapping/2D FFT calculations
# to occur at once. Make sure it's not on to begin with, or we'll never start.
LOCKDIR=$pipeline_DIR/LOCKDIR
rm -rf $LOCKDIR


# Make error report file if not already in existence
if [ ! -f "$pipeline_DIR/Error_Reports/" ]; then
    mkdir -p $pipeline_DIR/Error_Reports/
    touch $pipeline_DIR/Error_Reports/Error_Report.txt
fi


if [ "$RUN" == "Sims_Run" ]; then
	
	# Assemble the LOS into one array over which you loop
	source ShowSumClass/AssembleLOS.sh $DIRname $sqdeg $z $cosmol $los_start $los_end

	echo "You have selected: 
		  Mock run: 		$sqdeg Sqdeg
		  Shape Noise: 		$SN
		  Masking: 			$mask
 		  Gal. density: 	$gpam gal/sq. arcmin
		  F-Scale:			$SS Pxls
		  Clip threshold:	$sigma sigma
		  Mock cosmology:   $cosmol
		  Low-redshift lim: $zlo
		  Hi-redshift lim:  $zhi
		  No. Athena Bins:  $ThBins
		  Output filename:	$name
		  Output Subdir:	$DIRname" 	

	
	final_dir=$pipeline_DIR/Correlation_Function/$DIRname/ThBins$ThBins
	timerfile=$final_dir/Timings.$name.$gpam"GpAM".$sigma"sigma".LOS$los_start"_"LOS$los_end.txt
	clipfilepath=$pipeline_DIR/Clipping_K/Clip_Thresholds/ClipThreshold_"$sigma"sigma."$sqdeg"Sqdeg.SS$SS


else
    # NEED TO EDIT THIS PART TO SET UP MULTIPLE ZBcut DIRECTORIES

	# Assemble the Fields into array over which you loop
	Loop_Array=()
	timerfile_subname="" # assemble from input fields	
	count_Fields=0
	echo "You have selected Field(s):" 
	for f in ${Field[*]}; do
		Loop_Array+=("$f")
		echo "		$f"

		# Only record the G* part of the name, with no repeats.
		strip_name=$(echo ${f%_NOISE*})
		#echo "strip_name is "$strip_name
		#echo "timer_subname is "$timerfile_subname
		if [[ "$timerfile_subname" != *"$strip_name"* ]]; then
			timerfile_subname+="$strip_name"
			count_Fields=$((count_Fields+1))
		fi

	done

	NOISE_NLOS=0
	if [[ $DIRname == *"NOISE"* ]]; then 
		NOISE_NLOS=$((${#Field[@]} / count_Fields))
		timerfile_subname+="_NOISE"$NOISE_NLOS
	fi
	#echo $timerfile_subname

	echo "	Blind:		 	 $Blind"
	echo "	F-Scale:	 	 $SS Pxls,"
	echo "	Clip threshold:	 $sigma sigma"
	echo "  zlow             $zlo"
	echo "  zhigh:           $zhigh"
	echo "	Pure Noise?:		 $NOISE"
	echo "	No. of Noise Realisations: 		$NOISE_NLOS"


	final_dir=$pipeline_DIR/Correlation_Function/$DIRname/ThBins$ThBins
	timerfile=$final_dir/Timings.Blind$Blind.SS$SS.$sigma"sigma".$timerfile_subname.txt
	clipfilepath=$pipeline_DIR/Clipping_K/Clip_Thresholds/ClipThreshold_"$sigma"sigma.Blind$Blind.SS$SS

fi


# Make directories if they don't exist already                                                     
for DIR in $pipeline_DIR $data_DIR; do
    for dir in $DIRname $DIRname1 $DIRname2; do
        # DIRname1/2 only exist if working with multiple paramfiles/zbins                                   
        for subdir in "Mass_Recon" "Clipping_K" "Correlation_Function" "Tree_Correlation_Function"; do

	    make_directory="$DIR/$subdir/$dir"
	    if [[ "$make_directory" == *"Correlation_Function"* ]]; then
		make_directory=$make_directory/ThBins$ThBins
	    fi
	    	
	    if [ ! -d "$make_directory" ]; then
                mkdir -p $make_directory
            fi

        done
    done
done




echo "Master will loop over the following Fields/Lines of sight"
for f in ${Loop_Array[*]}; do
	echo $f
done


# rsync all required data to data_DIR if pipeline_DIR != data_DIR
if [ "$pipeline_DIR" != "$data_DIR" ]; then	source $pipeline_DIR/ShowSumClass/Change_DIR_of_outputs.sh; fi

# If an extra parameter file is submitted, & we are working with SLICS,
# and we are working with K1000 mocks, we need to run the above script twice.
# This is because this script copies SLICS cat's to the workers, and in the case of
# K1000 mocks, there's different cats for different zbins.
# Hence two different sets of cats need to be copied.
if [[ "$5" == *"param_files"* ]]; then
    echo "zlo2 and zhi2 are $zlo2 and $zhi2"
    zlo=$zlo2 # redefine these as those corresponding to 2nd paramfile (get set in Assemble_Combine-zbin-DIRname.sh)
    zhi=$zhi2 # (these are used by Change_DIR_...sh to find the correct cats to copy.
    source $pipeline_DIR/ShowSumClass/Change_DIR_of_outputs.sh
fi


rm -f $timerfile # remove the one already there, don't just peg onto old one.

echo "Time at which pipeline started: $(date)" >> $timerfile

# Now loop over the Fields/LOS
for Part in I II; do  # I II 
	counter=0
	for f in ${Loop_Array[*]}; 
	do


		# Wait if running n processes at once.
		if [ "$counter" == 1 ]; then 
			echo "Waiting for cores to clear."
			wait 
			counter=0
		fi	



		echo ""
		echo "Running pipeline part $Part on Field $f."
		echo "Final Field is ${Loop_Array[-1]}." 
		echo ""

		echo "Running pipeline part $Part on Field $f at $(date)" >> $timerfile




		if [ "$RUN" == "Sims_Run" ]; then
		        
			$pipeline_DIR/SkelePipeline_Part"$Part".sh $RUN $paramfile $f $los_end $paramfile2 &
			
			# Exit if there's a problem - NOTE. Doesnt always work, since it may take time to hit an error
			# ...whereas the pipeline has already moved onto commence subsequent LOS runs.
			if [ $? -eq 1 ]; then
				echo "$(date): SkelePipeline_Part$Part.sh failed for Field $f. Exiting Master."
				printf "$(date): SkelePipeline_Part$Part.sh failed for Field $f. Exiting Master." > $pipeline_DIR/Error_Reports/Error_Report.txt
				printf " -----------------------------------------------------------------------" > $pipeline_DIR/Error_Reports/Error_Report.txt
				wait # for other LOS to finish/fail
	    		exit 1
			fi
		else
			$pipeline_DIR/SkelePipeline_Part"$Part".sh $RUN $paramfile $f &
			if [ $? -eq 1 ]; then
				echo "$(date): SkelePipeline_Part$Part.sh failed for Field $f. Exiting Master."
				printf "$(date): SkelePipeline_Part$Part.sh failed for Field $f. Exiting Master." > $pipeline_DIR/Error_Reports/Error_Report.txt
				printf " -----------------------------------------------------------------------" > $pipeline_DIR/Error_Reports/Error_Report.txt
				wait # for other LOS to finish/fail
	    		exit 1
			fi
		fi




		counter=$((counter+1));
	done

	echo "Waiting to begin Part $Part of the pipeline..."
	wait
done




wait




echo "Time at which pipeline ended: $(date)" >> $timerfile



# 8. Plot the CF averaged across the fields. Give it all the inputs you gave this code.
echo "Executing plot_CorrFun.py"
#python $pipeline_DIR/Correlation_Function/plot_CorrFun.py $1 $2 $3 $4 $5 $6 $7 


if [ $? -eq 1 ]; then
    echo 'plot_CorrFun.py failed. Exiting Master.'
    exit 1
else
 	echo 'plotting ran correctly: Correlation plots and Covariance matrices made! (Unless its commented out).'
fi



# If results are being saved to data_DIR, rsync them back to pipeline_DIR
if [ "$pipeline_DIR" != "$data_DIR" ]; then     source $pipeline_DIR/ShowSumClass/Copy_Results_2PipelineDIR.sh; fi


if [ "$RUN" == "Sims_Run" ] && [[ $data_DIR == "/data/"* ]]; then
    #rm -rf $data_DIR/Mass_Recon/$DIRname/*_Xm_Ym_e1_e2_w.dat
    #rm -f $data_DIR/Mass_Recon/$DIRname/*SS${SS}*Ekappa.fits
    echo "Not deleting catalogue in Mass_Recon this time... do this by hand"
fi



# Change the name of the clip threshold so doesn't get used again
if [ "$sigma" != "X1" ] && [ "$sigma" != "X2" ] && [ "$sigma" != "X3" ]; then
	mv $clipfilepath.txt $clipfilepath.${Loop_Array[0]}to${Loop_Array[-1]}.txt
fi	






