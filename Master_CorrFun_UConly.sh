#!/bin/bash
# 05/11/15, B. M. Giblin, PhD student Edinburgh

# The Master Script loops over no. of LOS/Fields running 
# SkelePipeline.sh for each one <-- that contains the breakdown of 
# the clipping pipeline. 



# check validity of command line arguments
# Allow for max no. of args: KiDS_Run, all fields, Blind, SS, sigma
source ShowSumClass/FilterInputArgs.sh $1 $2 $3 $4 $5 $6 $7 $8 $9 
overall_DIR=$PWD



# The locking mechanism to prevent multiple mass mapping/2D FFT calculations
# to occur at once. Make sure it's not on to begin with, or we'll never start.
LOCKDIR=$overall_DIR/LOCKDIR
rm -rf $LOCKDIR


# Make error report file if not already in existence
if [ ! -f "$overall_DIR/Error_Reports/Error_Report.txt" ]; then
	mkdir -p $overall_DIR/Error_Reports/Error_Report.txt
fi


if [ "$RUN" == "Sims_Run" ]; then
	
	# Assemble the LOS into one array over which you loop
	source ShowSumClass/AssembleLOS.sh $DIRname $sqdeg $z $cosmol $los_start $los_end

	
	# Make Field subdirectories if they don't exist already
	if [ ! -d "$overall_DIR/Mass_Recon/$DIRname" ]; then
		mkdir $overall_DIR/Mass_Recon/$DIRname
	fi


	if [ ! -d "$overall_DIR/Clipping_K/$DIRname" ]; then
		mkdir $overall_DIR/Clipping_K/$DIRname
	fi


	if [ ! -d $overall_DIR/Correlation_Function/$DIRname ]; then
		mkdir $overall_DIR/Correlation_Function/$DIRname
	fi
	
	if [ ! -d $overall_DIR/Correlation_Function/$DIRname/ThBins$ThBins ]; then
		mkdir $overall_DIR/Correlation_Function/$DIRname/ThBins$ThBins
	fi

	if [ ! -d $overall_DIR/Tree_Correlation_Function/$DIRname/ThBins$ThBins ]; then
		mkdir -p $overall_DIR/Tree_Correlation_Function/$DIRname/ThBins$ThBins
	fi




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

	
	final_dir=$overall_DIR/Correlation_Function/$DIRname/ThBins$ThBins
	timerfile=$final_dir/Timings.$name.$gpam"GpAM".$sigma"sigma".LOS$los_start"_"LOS$los_end.txt
	clipfilepath=$overall_DIR/Clipping_K/Clip_Thresholds/ClipThreshold_"$sigma"sigma."$sqdeg"Sqdeg.SS$SS


else

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


	# Make Field subdirectories if they don't exist already
	if [ ! -d "$overall_DIR/Mass_Recon/$DIRname" ]; then
		mkdir $overall_DIR/Mass_Recon/$DIRname
	fi


	if [ ! -d "$overall_DIR/Clipping_K/$DIRname" ]; then
		mkdir $overall_DIR/Clipping_K/$DIRname
	fi


	if [ ! -d $overall_DIR/Correlation_Function/$DIRname ]; then
		mkdir $overall_DIR/Correlation_Function/$DIRname
	fi

	if [ ! -d $overall_DIR/Correlation_Function/$DIRname/ThBins$ThBins ]; then
		mkdir $overall_DIR/Correlation_Function/$DIRname/ThBins$ThBins
	fi

	if [ ! -d $overall_DIR/Tree_Correlation_Function/$DIRname/ThBins$ThBins ]; then
		mkdir -p $overall_DIR/Tree_Correlation_Function/$DIRname/ThBins$ThBins
	fi


	final_dir=$overall_DIR/Correlation_Function/$DIRname/ThBins$ThBins
	timerfile=$final_dir/Timings.Blind$Blind.SS$SS.$sigma"sigma".$timerfile_subname.txt
	clipfilepath=$overall_DIR/Clipping_K/Clip_Thresholds/ClipThreshold_"$sigma"sigma.Blind$Blind.SS$SS

fi


echo "Master will loop over the following Fields/Lines of sight"
for f in ${Loop_Array[*]}; do
	echo $f
done






rm -f $timerfile # remove the one already there, don't just peg onto old one.





echo "Time at which pipeline started: $(date)" >> $timerfile

# Now loop over the Fields/LOS
counter=0
for f in ${Loop_Array[*]}; 
do

	# Wait if running n processes at once.
	if [ "$counter" == 5 ]; then 
		echo "Waiting for cores to clear."
		wait 
		counter=0
	fi	



	echo ""
	echo "Running pipeline on Field $f."
	echo "Final Field is ${Loop_Array[-1]}." 
	echo ""

	echo "Running pipeline on Field $f at $(date)" >> $timerfile




	if [ "$RUN" == "Sims_Run" ]; then
		./SkelePipeline_UConly.sh $RUN $paramfile $f $los_end &
			
		# Exit if there's a problem - NOTE. Doesnt always work, since it may take time to hit an error
		# ...whereas the pipeline has already moved onto commence subsequent LOS runs.
		if [ $? -eq 1 ]; then
			echo "$(date): SkelePipeline_UConly.sh failed for Field $f. Exiting Master."
			printf "$(date): SkelePipeline_UConly.sh failed for Field $f. Exiting Master." > $overall_DIR/Error_Reports/Error_Report.txt
			printf " -----------------------------------------------------------------------" > $overall_DIR/Error_Reports/Error_Report.txt
			wait # for other LOS to finish/fail
    		exit 1
		fi
	else
		./SkelePipeline_UConly.sh $RUN $paramfile $f &
		if [ $? -eq 1 ]; then
			echo "$(date): SkelePipeline_UConly.sh failed for Field $f. Exiting Master."
			printf "$(date): SkelePipeline_UConly.sh failed for Field $f. Exiting Master." > $overall_DIR/Error_Reports/Error_Report.txt
			printf " -----------------------------------------------------------------------" > $overall_DIR/Error_Reports/Error_Report.txt
			wait # for other LOS to finish/fail
    		exit 1
		fi
	fi




	counter=$((counter+1));
done






wait




echo "Time at which pipeline ended: $(date)" >> $timerfile


if [ "$RUN" == "Sims_Run" ]; then
	#rm -f $overall_DIR/Mass_Recon/$DIRname/*.dat
	echo "Not deleting catalogue in Mass_Recon and Correlation_Function: handled by SkelePipeline_UConly.sh"
fi



# Change the name of the clip threshold so doesn't get used again
if [ "$sigma" != "X1" ] && [ "$sigma" != "X2" ] && [ "$sigma" != "X3" ]; then
	mv $clipfilepath.txt $clipfilepath.${Loop_Array[0]}to${Loop_Array[-1]}.txt
fi	








