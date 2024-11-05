#!/bin/bash

# 12/09/2016, B. M. Giblin, PhD student Edinburgh
# Filter the input given to bash scripts in the clipping pipeline

RUN=$1
if [ $# -eq 4 ] && [ "$1" == "Sims_Run" ]; then
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
	elif [ "$cosmol" == "fid" ]; then
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
	
	if [[ $Sys == *"IA"*  ||  $Sys == *"Bary"* ]]; then Prepend=${Prepend}${Sys}_; fi

	name=$name_start$name_end
	DIRname="$Prepend$sqdeg"Sqdeg$DIRname_start$DIRname_mid"$gpam"GpAM_z"$z"_ZBcut"$ZBcut$DIRname_end"
	z=$z$DIRname_end
	#echo $name
	#echo "DIRname is $DIRname"





elif [ $# -gt 2 ] && [ "$1" == "KiDS_Run" ]; then

	# See if there are INTEGERS at the end of Field arguments
	# If so, these are the PURE NOISE REALISATION Numbers (used to correct for KiDS-450 mask)
	if [[ "${@: -1}" =~ ^-?[0-9]+$ ]]; then
		los_end="${@: -1}"
	else
		los_end="-1"
	fi


	if [[ "${@:(-2):1}" =~ ^-?[0-9]+$ ]]; then
		los="${@:(-2):1}"
	else
		los="-1"
	fi

	Noise_Realisation_Numbers=($(seq $los 1 $los_end))

	paramfile=$2
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


	Blind=${args[0]}
	SS=${args[1]} # The no. of gal per square arcmin
	sigma=${args[2]}
	zlo=${args[3]}
	zhi=${args[4]}
	ThBins=${args[5]}
	MRres=${args[6]} # resolution to use in mass reconstruction (mask with this res should exist)
	OATH=${args[7]}
	NOISE=${args[8]}
	mask_variable=${args[9]}



	# Identify which of the inputs are Fields
	Field=()
	for f in $3 $4 $5 $6 $7 $8 $9; do
		if [ ${f:0:1} == "G" ]; then

			# Is this a Pure-Noise KiDS-run?
			if [ $los -ge 0 ] && [ $los_end -ge 0 ] && [ $NOISE == "Y" ]; then
				for nn in ${Noise_Realisation_Numbers[*]}; do Field+=($f"_NOISE"$nn); done
			else
				Field+=("$f")
			fi

		fi
	done

	if [ ${#Field[@]} -eq 0 ]; then 
		echo "No valid KiDS fields found"
		exit 1
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
	if [ "$MRres" == "" ] || [ "$MRres" == "-" ] || [ "$MRres" == "1arcmin" ] || [ "$MRres" == "arcmin" ]; then
		Prepend="" # Use fiducial value of 1arcmin for KiDS Run, do not append DIRname
		MRres="arcmin" # Set to this so mask can be easily identified.
	else
		Prepend="MRres"$MRres"_" # Not the default mass recon resolution. Prepend this to DIRname
	fi 


	# Check if you're running on shape noise ONLY
	if [ "$NOISE" == "N" ] || [ "$NOISE" == "" ] || [ "$NOISE" == "-" ]; then
		Add_Noise="" # Use fiducial value of 1arcmin for KiDS Run, do not append DIRname
		NOISE="N"
	else
		Add_Noise="_NOISE" # Not the default mass recon resolution. Prepend this to DIRname
		NOISE="Y"
	fi	


	# Check if we are masking this KiDS_Run
	if [ "$mask_variable" == "mask" ] || [ "$mask_variable" == "" ] || [ "$mask_variable" == "-" ]; then
		mask_variable="" # Use fiducial value of 1arcmin for KiDS Run, do not append DIRname
		Add_Mask=""
	else
		mask_variable="-nomask"
		Add_Mask="_NoMask" # Not the default mass recon resolution. Prepend this to DIRname
	fi		


	DIRname="KiDS_Fields_ZBcut$ZBcut$Add_Noise$Add_Mask"
	echo "DIRname is $DIRname"
	#echo $ThBins



else
  echo "ARGUMENTS:"
  echo "	1. Sims_Run"
  echo "	1. <address of param file>"
  echo "	2. current los"
  echo "	3. end LOS"
  echo "	----------------OR----------------"
  echo "	1. KiDS_Run"
  echo "	2. <address of param file>"
  echo "	3-7. etc. G.., G.., G.., (KiDS Field names, as many or as little as desired)"
  echo "    8. Noise los number"
  echo "    9. Noise los_end number"
  echo "	INVALID INPUT. EXITING CODE"
  exit 1
fi



