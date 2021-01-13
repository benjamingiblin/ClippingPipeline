#!/bin/bash
# 03/11/15, B. M. Giblin, PhD student Edinburgh

# Runs Athena, calculating the shear-shear (i.e. cosmological lensing) correlation function, as a function on theta (arcmin). 
# It relies on the config_tree configuration file. Outputs a catalogue.
# Input for this code is a catalogue of format: X(deg) Y(deg) gamma1 gamma2 weight

# Must be ran from the directory above.



source ShowSumClass/FilterInputArgs.sh $1 $2 $3 $4 $5
overall_DIR=$PWD
direct=$overall_DIR/Correlation_Function
athena=/home/bengib/athena_1.7/src/athena



if [ "$RUN" == "Sims_Run" ]; then
	start_name=$direct/$DIRname/$name."$gpam"GpAM.LOS$los
	start_name_out=$direct/$DIRname/ThBins$ThBins/$name."$gpam"GpAM.LOS$los
	
	cn=$(echo $(($los%8)) | bc -l) 
	# cn distinguishes config files used for simultaneous Athena runs
	# take remainder of LOS/8 
	output_message="$sqdeg Sqdeg LOS $los"
	if [ "$sqdeg" == "60" ]; then
		thmin=0.5 
		thmax=300. # in arcmin
		radec=0 # cartesian
	elif [ "$sqdeg" == "100" ]; then
		thmin=0.5 
		thmax=300. # in arcmin
		radec=0 #cartesian
	else
		thmin=0.5
		thmax=300. # in arcmin
		radec=1 #SPHERICAL
	fi

else
	start_name=$direct/$DIRname/$Field.Blind$Blind
	start_name_out=$direct/$DIRname/ThBins$ThBins/$Field.Blind$Blind
	cn=$Field	
	output_message='Field '$Field
	thmin=0.5
	thmax=300. # in arcmin
	radec=1 # spherical	 
fi



# define difference in filenames between clipped and unclipped runs
rc_middle=SS$SS.rCLIP_"$sigma"sigma
uc_middle=ORIG

for type in $uc_middle $rc_middle; do

	infile=$start_name.$type.ThetaX_ThetaY_e1_e2_w.Std.asc
	outfile=$start_name_out.$type.CorrFun.asc

	ShowSumClass/AssembleConfigTree.sh $cn $infile $infile $thmin $thmax $ThBins $radec $OATH

	if [ "$type" == $uc_middle ]; then #&& [ ! -f $start_name_out.ORIG.CorrFun.asc ]; then
									   # make it so it overwrites unclipped? comment out &&...
		echo "$(date): Running Athena on Unclipped for "$output_message
		$athena -q -c $direct/config_files/config_tree$cn --out_xi $outfile


	elif [ "$type" == $rc_middle ]; then

		echo "$(date): Running Athena on clipped for "$output_message
		$athena -q -c $direct/config_files/config_tree$cn --out_xi $outfile

	fi

	# Check if Athena failed. Save an error message if so.
	if [ $? -eq 1 ]; then
		echo "$(date): Athena failed when computing $outfile" > $direct/$DIRname/ThBins$ThBins/Error_Report.txt
	else
		echo "Athena was successful for "$output_message
	fi



done



# Calculate the cross-covariance between unclipped and delta-clip?
cross_corr='Y'
if [ "$cross_corr" == 'Y' ]; then
	
	echo "$(date): Running Athena for cross-correlation between unclipped and delta-clip for "$output_message
	infile1=$start_name.ORIG.ThetaX_ThetaY_e1_e2_w.Std.asc
	infile2=$start_name.SS$SS.deltaCLIP_"$sigma"sigma.ThetaX_ThetaY_e1_e2_w.Std.asc
	outfile=$start_name_out.SS$SS.deltaXOrigCLIP_"$sigma"sigma.CorrFun.asc
	ShowSumClass/AssembleConfigTree.sh $cn $infile1 $infile2 $thmin $thmax $ThBins $radec $OATH
	$athena -q -c $direct/config_files/config_tree$cn --out_xi $outfile


	# Check if Athena failed. Save an error message if so.
	if [ $? -eq 1 ]; then
		echo "$(date): Athena failed when computing $outfile" > $direct/$DIRname/ThBins$ThBins/Error_Report.txt
	else
		echo "Athena was successful for "$output_message
	fi


	echo "$(date): Running Athena for auto-correlation on delta-clip for "$output_message
	infile1=$start_name.SS$SS.deltaCLIP_"$sigma"sigma.ThetaX_ThetaY_e1_e2_w.Std.asc
	infile2=$start_name.SS$SS.deltaCLIP_"$sigma"sigma.ThetaX_ThetaY_e1_e2_w.Std.asc
	outfile=$start_name_out.SS$SS.deltaCLIP_"$sigma"sigma.CorrFun.asc
	ShowSumClass/AssembleConfigTree.sh $cn $infile1 $infile2 $thmin $thmax $ThBins $radec $OATH
	$athena -q -c $direct/config_files/config_tree$cn --out_xi $outfile

	# Check if Athena failed. Save an error message if so.
	if [ $? -eq 1 ]; then
		echo "$(date): Athena failed when computing $outfile" > $direct/$DIRname/ThBins$ThBins/Error_Report.txt
	else
		echo "Athena was successful for "$output_message
	fi

fi



#rm -f $direct/config_files/config_tree$cn



# remove some of the annoying things athena creates
rm -f $overall_DIR/xi.resample*
rm -f $overall_DIR/log
rm -f $direct/xi.resample*
rm -f $direct/log






