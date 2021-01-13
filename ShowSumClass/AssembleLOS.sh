#!/bin/bash

# 03/11/16, B. M. Giblin, PhD Student Edinburgh
# This code is ran by the Master -> It assembles the mocks, extracts the LOS no.s
# And returns a list to the Master.


data_DIR='/data/bengib/Clipping_SimsLvW/'
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
	    
	    if [ "$z" == "KiDS1000" ]; then
		echo "THIS IS A SLICS KiDS1000-like RUN..."
		
		# Annoying, the redshift cuts have already been made on KiDS1000 mocks, so
		# we're limited to using catatalogues with the KiDS1000 bins:
		if [ "$zlo" == "0.1" ] && [ "$zhi" == "0.3" ]; then
                    bin_name="bin1"
		elif [ "$zlo" == "0.3" ] && [ "$zhi" == "0.5" ]; then
                    bin_name="bin2"
		elif [ "$zlo" == "0.5" ] && [ "$zhi" == "0.7" ]; then
                    bin_name="bin3"
		elif [ "$zlo" == "0.7" ] && [ "$zhi" == "0.9" ]; then
                    bin_name="bin4"
		elif [ "$zlo" == "0.9" ] && [ "$zhi" == "1.2" ]; then
                    bin_name="bin5"
		else
                    echo "You have set z to KiDS1000. In this case, zlo & zhi must be consecutive numbers from:"
                    echo " 0.1, 0.3, 0.5, 0.7, 0.9, 1.2 (i.e. the bins used in KiDS1000 cosmic shear)."
                    echo "But you have set zlo and zhi to $zlo, $zhi. Not programmed to deal with this shit. EXITING."
                    exit
		fi
		filename=/disk10/jharno/MockProducts/KiDS1000/${bin_name}/GalCatalog_KiDS1000_${bin_name}_LOS
	    else
		echo "THIS IS A SLICS KV450-like RUN..."
		filename=/disk10/jharno/MockProducts/KV450/GalCatalog_KV450_LOS
	    fi
	    
	    
	elif [[ "$z" == *"LSST"* ]]; then
	    echo "THIS IS A SLICS Generic (LSST-like) RUN..."
	    filename=/disk10/jharno/MockProducts/Generic/GalCatalog_LOS
	fi
	
	Files=${filename}*.fits
        # Find the LOS that exist
	lines_of_sight=()
	for f in $Files; do
	    f2=$(echo $f | rev | cut -d_ -f1 | rev) # gives LOS***.fits
	    f3=$(echo $f2 | rev | cut -d. -f2 | rev) # gives LOS***
	    f4=$(echo $f3 | rev | cut -dS -f1 | rev) # gives ****
	    lines_of_sight+=("$f4")
	done
    fi

    if [[ $DIRname == *"SNCycle"* ]]; then
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





# Only try to assemble runs into order if you're NOT doing a noise-cycle run.
if [[ $DIRname != *"SNCycle"* ]]; then
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


