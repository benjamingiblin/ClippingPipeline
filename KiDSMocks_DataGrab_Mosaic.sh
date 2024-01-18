#!/bin/bash
# 15/12/2022, B. Giblin, Postdoc Barcelona
# The format of the KiDS-1000 mosaic mocks is so different
# to all of the other (cosmo)-SLICS mocks that we have this
# separate script just to access them and pull out the rel. columns

pipeline_DIR='/home/bengib/Clipping_Pipeline/'
data_DIR='/data/bengib/Clipping_Pipeline/'

source $pipeline_DIR/ShowSumClass/FilterInputArgs.sh $1 $2 $3 $4 $5

# Get the tomo bin
source $pipeline_DIR/ShowSumClass/Identify_KiDS1000_zbin.sh $zlo $zhi
bnum=${bin_name#"bin"}
# Get the region R
R=${los#*R}
R=${R%n*}
# And the LOS...
los_num=${los%R*}

# Get the end part of DIRname
IFS='_' read -ra ENDname <<< "$DIRname"
if [[ ${ENDname[-1]} == *"Cosmol"* ]]; then
    # It's a cosmoSLICS mosaic run...
    
    # Get the cosmol ID by splitting on 'Cosmol'
    IFS="Cosmol" read -ra ID <<< "${ENDname[-1]}"
    cosmol=${ID[-1]}  # Only want the last element.                                                                                            
    cosmol_fname=$cosmol
    if [ "$cosmol_fname" != "fid" ]; then # Then it's an integer, make it 2sigfig
	while [[ ${#cosmol_fname} -lt 2 ]] ; do cosmol_fname="0${cosmol_fname}"; done
    fi

    # If los<26 read in 'a' seeds, otherwise 'f' seeds
    if [ "$los_num" -lt "26" ]; then
        seed="a"
	los_col=$los_num
    else
        seed="f"
	los_col=$((los_num-25))
    fi

    mosaic_DIR=/disk10/jharno/PeakCount/Mocks/Cosmo_KiDS1000/${cosmol_fname}_${seed}/
    input=${mosaic_DIR}/AllLOS/KiDS1000_MocksCat_${cosmol_fname}_${seed}_${bnum}_LOS1-25_R${R}.dat

    # Use the LOS to determinewhich columns to pull out:
    col_g1=$((los_col-1))  # subtract 1 because cosmoSLICS LOS start at 1, not 0.
    col_g1=$((col_g1*3))   # x3 because each LOS has (g1,g2,k) saved in cat
    col_g1=$((col_g1+9))   # add 8 as there's 8 (ra,dec,etc...) columns at start

else
    # It's a SLICS mosaic KiDS1000 run...
    # multiple LOS stored in a single cat, find which to copy:
    if [ "$los_num" -lt "74" ]; then
	echo "EXITING BECAUSE LOS $los_num IS LESS THAN LOS_MIN 74"
	exit
    elif [ "$los_num" -lt "200" ]; then
        los_low=74
    elif [ "$los_num" -lt "299" ]; then
        los_low=200
    fi
    mosaic_DIR=/disk10/jharno/PeakCount/Mocks/Cov_KiDS1000/SLICS
    input=${mosaic_DIR}/KiDS1000_MocksCat_SLICS_${bnum}_LOS${los_low}*_R${R}.dat

    # Use the LOS to determinewhich columns to pull out:
    col_g1=$((los_num-los_low))
    col_g1=$((col_g1*3))
    col_g1=$((col_g1+9))
    
fi

col_g2=$((col_g1+1))
echo "col_g1 is $col_g1, col_g2 is $col_g2"

output=$data_DIR/Mass_Recon/$DIRname/$name."$gpam"GpAM.LOS"$los"_Xm_Ym_e1_e2_w.dat
if [ ! -d "$data_DIR/Mass_Recon/$DIRname" ]; then
    # This will only be the case if we're doing KiDS-1000 run no-tomography
    # in which we do an impromptu scroll through 5 tomo bins.
    mkdir -p $data_DIR/Mass_Recon/$DIRname
fi

# pull out RA, Dec, g1, g2, w, z, m_Angus
# (weight is not unity with mosaic mocks)
awk -v cg1=$col_g1 -v cg2=$col_g2 '{print $1, $2, $cg1, $cg2, $5, $8}' < $input > $output

