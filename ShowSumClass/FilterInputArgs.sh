#!/bin/bash

# 12/09/2016, B. M. Giblin, PhD student Edinburgh
# Filter the input given to bash scripts in the clipping pipeline

RUN=$1
if [ "$RUN" == "Sims_Run" ] || [ "$RUN" == "KiDS_Run" ]; then
	paramfile=$2
	los=$3
	los_start=$3 # Variable used only by Master
	los_end=$4

	c=0
	args=()
	while read p
 	 do
 		# $p is the first line of paramfile
		IFS='#' read -r id string <<< "$p" 		# split $p on '#' delim and put
												# first elem into $id
		id="$(echo -e "${id}" | tr -d '[[:space:]]')" # strip the white space from string

		if [ "$id" != "" ]; then
			#echo $id 
    		args[c]=$(echo "$id")
    		c=$((c+1));
		fi
	done < $paramfile


	gpam=${args[0]}
	SS=${args[1]} # The no. of gal per square arcmin
	sigma=${args[2]}
	SN=${args[3]}
	mask=${args[4]}
	z=${args[5]}
	PS=${args[6]}
	sqdeg=${args[7]}
	cosmol=${args[8]}
	zlo=${args[9]}
	zhi=${args[10]}
	ThBins=${args[11]}
	MRres=${args[12]} # resolution to use in mass reconstruction (mask with this res should exist)
	OATH=${args[13]}
	Sys=${args[14]}
	#echo "Sys is $Sys"
	#echo $gpam, $SS, $sigma, $SN, $mask, $z, $PS, $sqdeg, $cosmol, "END"


	# Now assemble the name and DIRname

	# SHAPE NOISE	
	if [ $SN == "ALL" ] || [ $SN == "All" ] || [ $SN == "all" ]; then
		name_start='NOISE_'
		DIRname_start='_NOISE'
	elif [ $SN == "ALL-KiDS" ]; then
		name_start='NOISE-KiDS_'
		DIRname_start='_NOISE-KiDS'
	elif [[ $SN == *"Cycle"* ]] || [[ $SN == *"cycle"* ]]; then

	    if [[ "$z" == *"KiDS1000"* ]]; then
		# Use KiDS1000 redshift binning and hence SN values.
		# These are not saved in the paramfile when running SNCycle,
		# Use following code to identify them:
		source $pipeline_DIR/ShowSumClass/Identify_KiDS1000_zbin.sh $zlo $zhi
		name_start="SN${sigma_e}_"
		
	    else 
		name_start='SN0.28_'
	    fi
	    DIRname_start='_SNCycle'
	    
	elif [ $SN == "0." ] || [ $SN == "0" ]; then
		name_start='NF_'
		DIRname_start='_NF'
	else
		name_start='SN'$SN'_'
		DIRname_start='_SN'$SN
	fi



	
	# MASKING
	if [ $mask == "nomask" ]; then
		mask_variable='-nomask'
		name_end='test'
		DIRname_mid='_NoMask_'

	elif [ $mask == "mask" ]; then
		mask_variable=''
		name_end='Mask'
		DIRname_mid='_Mask_'

	elif [ $mask == "G9mask" ]; then
		mask_variable=''
		name_end='G9Mask'
		DIRname_mid='_G9Mask_'

	elif [ $mask == "W3mask" ]; then
		mask_variable=''
		name_end='W3Mask'
		DIRname_mid='_W3Mask_'
		
	elif [[ $mask == *"mosaic"* ]]; then
	    mask_variable=''
	    name_end='Mosaic'
	    DIRname_mid='_Mosaic_'
	else
		echo "mask variable in paramfile must be 'mask' or 'nomask'. Please fix this."
		exit 1
	fi



	# Check if high/low sigma8
	if [ "$cosmol" == "high_sigma8" ]; then
		DIRname_end='_HS8'
	elif [ "$cosmol" == "low_sigma8" ]; then
	    DIRname_end='_LS8'
	elif [ "$cosmol" == "fid" ] || [ "$cosmol" == "KiDS1000" ]; then
	    DIRname_end="_Cosmol$cosmol"
	elif [ "$cosmol" -eq "$cosmol" ] 2>/dev/null; then
		# If cosmol is an integer, you're either running DH10 or cosmoSLICS mocks		
		DIRname_end="_Cosmol$cosmol"
	else
		DIRname_end=''
	fi



	# check if there's a zB cut
	# first check zlo and zhi are valid numbers
	zlo_variable=$(echo $zlo | grep -Eq '^[-+]?[0-9]+\.?[0-9]*$' && echo Match)
	zhi_variable=$(echo $zhi | grep -Eq '^[-+]?[0-9]+\.?[0-9]*$' && echo Match)
	if [ "$zlo_variable" == "Match" ] && [ "$zhi_variable" == "Match" ] && [ $(echo "$zhi > $zlo" | bc) -eq 1 ]; then
		ZBcut="$zlo-$zhi"
	else
		ZBcut="None"
	fi


	# check number of Athena theta bins is an integer
	if [ "$ThBins" -eq "$ThBins" ] 2>/dev/null; then
		echo ""
	else
		echo "The number of Athena theta bins in param file is not an integer. \n Setting this number to default of 9 bins"
		ThBins=9
	fi
		

	# check what resolution you want to use in mass reconstruction
	if [ "$sqdeg" == "100" ] || [ "$sqdeg" == "60" ] || [ "$sqdeg" == "36" ]; then
		if [ "$MRres" == "" ] || [ "$MRres" == "-" ] || [ "$MRres" == "5arcs" ]; then
			Prepend="" # Use fiducial value of 5arcs, do not append to DIRname
		else
			Prepend="MRres"$MRres"_" # Not the default mass recon resolution. Prepend this to DIRname
		fi
	elif  [ "$sqdeg" == "5000" ]; then
		Prepend="NSIDE"$MRres"_"
	fi
	
	if [[ $Sys == *"IA"*  ||  $Sys == *"Bary"* || $Sys == *"dz"* || $Sys == *"SLC"* ]]; then Prepend=${Prepend}${Sys}_; fi

	name=$name_start$name_end
	DIRname="$Prepend$sqdeg"Sqdeg$DIRname_start$DIRname_mid"$gpam"GpAM_z"$z"_ZBcut"$ZBcut$DIRname_end"
	###if [ "$RUN" == "KiDS_Run" ]; then DIRname=KiDS1000Data_${DIRname}; fi

	z=$z$DIRname_end
	#echo $name
	#echo "DIRname is $DIRname"


else
  echo "ARGUMENTS:"
  echo "	1. Sims_Run OR KiDS_Run"
  echo "	1. <address of param file>"
  echo "	2. current los"
  echo "	3. end LOS"
  echo "	INVALID INPUT. EXITING CODE"
  exit 1
fi



