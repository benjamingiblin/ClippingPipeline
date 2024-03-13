#!/bin/bash

# 03/11/16, B. M. Giblin, PhD Student Edinburgh
# This code is ran by the Master -> It assembles the mocks, extracts the LOS no.s
# And returns a list to the Master.


data_DIR='/data/bengib/Clipping_Pipeline/'
if [ ! -d "$data_DIR" ]; then mkdir $data_DIR; fi
	
DIRname=$1
sqdeg=$2
z=$3
cosmol=$4
LOS_start=$5
LOS_end=$6
echo "LOS...: $LOS_start $LOS_end"

if [ "$sqdeg" == "60" ]; then
	mocks60dir='/home/bengib/SLICS_60/'
	cd $mocks60dir/z_"$z"/
	Files="$z"gamma1_weight.dat_LOS*.fits

	# Find the LOS that exist
	lines_of_sight=()
	for f in $Files; do 
		f2=$(echo $f | rev | cut -d. -f2 | rev) # gives dat_LOS***
		f3=$(echo $f2 | rev | cut -d_ -f1 | rev) # gives LOS***
		f4=$(echo $f3 | rev | cut -dS -f1 | rev) # gives ***
		lines_of_sight+=("$f4")
	done







elif [ "$sqdeg" == "100" ]; then

    # Figure out if it's SLICS or cosmoSLICS using 'cosmol' variable
    # If cosmoSLICS, theres no missing LOS, and everything is simple:
    if [ "$cosmol" == "fid" ] || [ "$cosmol" -eq "$cosmol" ] 2>/dev/null; then
	echo "THIS IS A COSMO-SLICS RUN..."
	lines_of_sight=($(seq $LOS_start 1 $LOS_end))

    else
	# But if it's a SLICS run, there's missing LOS, so we have to find which ones exist:
	# First decide which SLICS: KiDS-like (KV450 or KiDS1000) or Generic (LSST-like)?
	if [[ "$z" == *"KiDS"* ]]; then
	    
	    if [[ "$z" == *"KiDS1000"* ]]; then
		echo "THIS IS A SLICS KiDS1000-like RUN..."
		# Annoying, the redshift cuts have already been made on KiDS1000 mocks, so
		# we have to access different directories for each zbin:
		source $pipeline_DIR/ShowSumClass/Identify_KiDS1000_zbin.sh $zlo $zhi
	     	filename=/disk10/jharno/MockProducts/KiDS1000/${bin_name}/GalCatalog_KiDS1000_${bin_name}_LOS
	    else
		echo "THIS IS A SLICS KV450-like RUN..."
		filename=/disk10/jharno/MockProducts/KV450/GalCatalog_KV450_LOS
	    fi
	    
	    
	elif [[ "$z" == *"LSST"* ]]; then
	    echo "THIS IS A SLICS Generic (LSST-like) RUN..."
	    filename=/disk10/jharno/MockProducts/Generic/GalCatalog_LOS
	fi

	# The SLICS mosaic mocks have different missing LOS - 7 missing in the 200-299 files
	# which are LOS 210-216 (mistake on JHD's part)
	# so manually adjust LOS_end to 292 if working with these mocks.
	# Means my LOS won't match JHD's, as mine are continuous. 
	# LOS 135 & 140 are missing between 74-199; so reduce LOS_end to 197 if in
	# this range to compensate.
	# MEANS THE MISSING LOS WILL BE (198,199) and (>292); i.e. the end of each 74-199 / 200-299 file.
	
	if [[ "$DIRname" == *"Mosaic"* ]]; then
	    if [ "$LOS_end" -gt "197" ] && [ "$LOS_end" -lt "200" ]; then
		echo " -x-x-x- MANUALLY SETTING LOS_end TO 197 TO AVOID ACCESSING MISSING LOS -x-x-x-"
		LOS_end=197
	    elif [ "$LOS_end" -gt "292" ]; then
		echo " -x-x-x- MANUALLY SETTING LOS_end TO 292 TO AVOID ACCESSING MISSING LOS -x-x-x-"
		LOS_end=292
	    fi

	    # Here we set up the LOS we'll use without using the non-mosiac file names.
	    # We're not using these because we're designing this so the missing LOS
	    # between the mosaic and non-mosaic mocks are different (this is the easiest way). 
	    lines_of_sight=()
	    for l in `seq $LOS_start $LOS_end`; do
		if [ "$l" -gt "197" ] && [ "$l" -lt "200" ]; then
		    echo "Skipping LOS $l" # skip (198,199)
		else
		    lines_of_sight+=("$l")
		fi
	    done
	    
	else
	    # non-mosaic run; use the filenames of SLICS to omit missing LOS.
	    Files=${filename}*.fits
            # Find the LOS that exist
	    lines_of_sight=()
	    for f in $Files; do
		f2=$(echo $f | rev | cut -d_ -f1 | rev) # gives LOS***.fits
		f3=$(echo $f2 | rev | cut -d. -f2 | rev) # gives LOS***
		f4=$(echo $f3 | rev | cut -dS -f1 | rev) # gives ****
		# following is only necessary for a mosaic (SLICS) run, but doesn't harm non-mosaic (SLICS) run
		if [ "$f4" -ge "$LOS_start" ] && [ "$f4" -le "$LOS_end" ]; then lines_of_sight+=("$f4"); fi
	    done
	    lines_of_sight=( $( printf "%s\n" "${lines_of_sight[@]}" | sort -n ) ) # sort into order
	fi
	
    fi


    if [[ $DIRname == *"Mosaic"* ]] && [[ $DIRname == *"SNCycle"* ]]; then
	# In this case we ened to cycle through both mosaic tiles and noise realisations...:
	
	tmp=${mask#*mosaic}         # returns e.g. mosaic tile IDs (X-Y)
	start_cycle_r=${tmp%-*}
	end_cycle_r=${tmp#*-}
	R_cycles=($(seq $start_cycle_r 1 $end_cycle_r))

	tmp=${SN#*ycle}         # returns e.g. 0-5
	start_cycle=${tmp%-*}
        end_cycle=${tmp#*-}
        noise_cycles=($(seq $start_cycle 1 $end_cycle))
	
	Loop_Array=()
	for l in ${lines_of_sight[*]}; do
	    for r in ${R_cycles[*]}; do
		for n in ${noise_cycles[*]}; do
		    Loop_Array+=("${l}R${r}n${n}")
		done
	    done
	done
    
    elif [[ $DIRname == *"SNCycle"* ]]; then
	tmp=${SN#*ycle} 	# returns e.g. 0-5
	start_cycle=${tmp%-*}
	end_cycle=${tmp#*-}
	noise_cycles=($(seq $start_cycle 1 $end_cycle))
	Loop_Array=()
	for l in ${lines_of_sight[*]}; do
	    for n in ${noise_cycles[*]}; do
		Loop_Array+=("${l}n${n}")
	    done
	done
	
    elif [[ $DIRname == *"Mosaic"* ]]; then
	tmp=${mask#*mosaic}         # returns e.g. mosaic tile IDs (X-Y)
	start_cycle_r=${tmp%-*}
	end_cycle_r=${tmp#*-}
        R_cycles=($(seq $start_cycle_r 1 $end_cycle_r))
	Loop_Array=()
        for l in ${lines_of_sight[*]}; do
            for r in ${R_cycles[*]}; do
                Loop_Array+=("${l}R${r}")
            done
        done
    fi
    



elif [ "$sqdeg" == "5000" ]; then

    lines_of_sight=($(seq $LOS_start 1 $LOS_end))



elif [ "$sqdeg" == "36" ]; then
    if [[ $DIRname == *"SNCycle"* ]]; then
	tmp=${SN#*ycle} 	# returns e.g. 0-5
	start_cycle=${tmp%-*}
	end_cycle=${tmp#*-}
	noise_cycles=($(seq $start_cycle 1 $end_cycle))
	los_cycles=($(seq $LOS_start 1 $LOS_end))
	Loop_Array=()
	for l in ${los_cycles[*]}; do
	    for n in ${noise_cycles[*]}; do
		Loop_Array+=("${l}n${n}")
	    done
	done
		

    else
	lines_of_sight=($(seq $LOS_start 1 $LOS_end))
    fi

else
	echo "		Pipeline currently only accepts 36 sqdeg, 60 sqdeg, 100 sqdeg and 5000 sqdeg mocks."
	echo "		You have (perhaps wrongly) specified $sqdeg sqdeg mocks." 
	echo "		If you want to persist, modify Sims_DataGrab.py & AssembleLOS.sh to accommodate this"
	exit 1
fi



# Only try to assemble runs into order if you're NOT doing a noise-cycle run OR mosaic run
if [[ $DIRname != *"SNCycle"* ]] && [[ $DIRname != *"Mosaic"* ]]; then
    # Sort the LOS array into ascending order:
    lines_of_sight=( $( printf "%s\n" "${lines_of_sight[@]}" | sort -n ) )

    # Of those that exist, select the ones in range: LOS_start --> LOS_end
    Loop_Array=()
    for l in ${lines_of_sight[*]}; do
	if [ "$l" -ge "$LOS_start" ] && [ "$l" -le "$LOS_end" ]; then Loop_Array+=("$l"); fi
    done


    cd $data_DIR

    if [ -f Mass_Recon/$DIRname/OkayLOS.txt ]; then
	rm -f Mass_Recon/$DIRname/OkayLOS.txt
    fi

    printf "%s\n" "${Loop_Array[@]}" > Mass_Recon/$DIRname/OkayLOS.txt


fi


