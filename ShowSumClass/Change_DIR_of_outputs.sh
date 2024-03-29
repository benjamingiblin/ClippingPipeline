#!/bin/bash

# 30/10/2017, B. M. Giblin, PhD Student, Edinburgh
# This code is useful if working on a worker of a supercomputer that is not directly connected to your pipeline_DIR.
# This code will only get executed if pipeline_DIR and data_DIR are not the same in Master_CorrFun_ByParts.sh
# This will copy the pipeline directory architecture, AND the data on which the pipeline runs
# FROM pipeline_DIR to data_DIR. 
# There is a LOCKDIR part to this so that you do not have lots of workers trying to copy data simulataneously,
# as this will only freeze up the connextors from workers - head - storage.

echo "You are working on ..."
hostname
module load compiler mkl			# make sure all necessary fortarn libraries are available.


# Firstly, make $data_DIR if it does not exist
if [ ! -d "$data_DIR" ]; then
	mkdir -p $data_DIR
fi

# Recreate the pipeline architecture but in data_DIR
if [ ! -d "$data_DIR/Mass_Recon/$DIRname" ]; then
	mkdir -p $data_DIR/Mass_Recon/$DIRname
fi


if [ ! -d "$data_DIR/Clipping_K/$DIRname" ]; then
	mkdir -p $data_DIR/Clipping_K/$DIRname
fi


if [ ! -d $data_DIR/Correlation_Function/$DIRname/ThBins$ThBins ]; then
	mkdir -p $data_DIR/Correlation_Function/$DIRname/ThBins$ThBins
fi
		

if [ ! -d $data_DIR/Tree_Correlation_Function/$DIRname/ThBins$ThBins ]; then
	mkdir -p $data_DIR/Tree_Correlation_Function/$DIRname/ThBins$ThBins
fi
		

# Where the data/masks/mocks are stored. These lines are automatically changed by Install_Pipeline.sh
KiDS_Data_DIR=/home/bengib/KiDS450/
KiDS_Mask_DIR=/home/bengib/KiDS450/
	
SLICS_KV450_DIR='/disk10/jharno/MockProducts/KV450/'
SLICS_KiDS1000_DIR='/disk10/jharno/MockProducts/KiDS1000/'
SLICS_Mosaic_DIR='/disk10/jharno/PeakCount/Mocks/Cov_KiDS1000/SLICS_Angus/'
SLICS_Generic_DIR='/disk10/jharno/MockProducts/Generic/'

cosmoSLICS_DIR=/home/jharno/Projects/cosmoSLICS/
cosmoSLICS_Mosaic_DIR=/disk10/jharno/PeakCount/Mocks/Cosmo_KiDS1000_Angus_m/
	
G9_Mask_DIR=/home/bengib/KiDS450/
W3_Mask_DIR=/home/bengib/WMAP_Masks/
	
MiraTitan_DIR=/home/bengib/MiraTitan/

DH10_DIR=/home/bengib/DH10_Mocks/FaLCoNS/

# Where the data/masks/mocks WILL BE stored. These lines are automatically changed by Install_Pipeline.sh
KiDS_Data_dataDIR=$data_DIR/KiDS450/
KiDS_Mask_dataDIR=$data_DIR/KiDS450/

SLICS_dataDIR='/data/bengib/Clipping_Pipeline//SLICS_100/'
cosmoSLICS_dataDIR='/data/bengib/Clipping_Pipeline//SLICS_100/'

G9_Mask_dataDIR=$data_DIR/KiDS450/
W3_Mask_dataDIR=$data_DIR/WMAP_Masks/

MiraTitan_dataDIR=$data_DIR/MiraTitan/

DH10_dataDIR=$data_DIR/DH10_Mocks/FaLCoNS/


# Now figure out which mocks/data/masks you are working with...
if [ "$RUN" == "Sims_Run" ]; then
	# DH10
	if [ "$sqdeg" == "36" ]; then
	    if [ ! -d "$DH10_dataDIR" ]; then mkdir -p $DH10_dataDIR; fi

	    # If dealing with the fiducial cosmology, only copy the necessary LOS directories
	    # (takes too long to copy all)
	    cosmol=$(echo ${DIRname##*Cosmol})
	    if [ $cosmol -eq 0 ]; then
		num_rt_cats=5
		for ll in `seq $los_start $los_end`; do
		    intlos=$(echo ${ll%n*})
		    realisation_number=$(echo "$((intlos / num_rt_cats + 1))")
		    # need to make realisation_number longer: 5 --> 005
		    while [[ ${#realisation_number} -lt 3 ]] ; do realisation_number="0${realisation_number}"; done
		    catalogue_number=$(echo "$((intlos % num_rt_cats))")
		    allcatalogues="${DH10_DIR}/Cosmol$cosmol/N*_L*_Om*_s8*_w-1.00-${realisation_number}/rt/gal/*/FaLCoNS_with_KiDS_nofz_0.5_z_0.9.cat"
		    for catcat in ${allcatalogues[*]}; do
			stripdirfile=$(echo ${catcat##*$DH10_DIR})                         # gets from /Cosmol0/... onwards to file
			stripdir=$(echo ${stripdirfile%/FaLCoNS_with_KiDS*})          # removes the file part
			if [ ! -d "$DH10_dataDIR/$stripdir" ]; then mkdir -p $DH10_dataDIR/$stripdir; fi
			echo "Copying $catcat to $DH10_dataDIR/$stripdir/"
			rsync -avz $catcat $DH10_dataDIR/$stripdir/
		    done
		done
	    else
		rsync -avz $DH10_DIR/Cosmol$cosmol $DH10_dataDIR
	    fi

	    if [ ! -d "$W3_Mask_dataDIR" ]; then mkdir -p $W3_Mask_dataDIR; fi
	    rsync -avz $W3_Mask_DIR/W3.16bit.5arcs.reg.Now_"$sqdeg"sqdeg.fits $W3_Mask_dataDIR/



	    
	# SLICS or cosmoSLICS
	elif [ "$sqdeg" == "100" ]; then
	    # Is it CosmoSLICS (Fid/intger) OR SLICS (Fiducial/high_sigma8/low_sigma8)?
	    IFS='_' read -ra ENDname <<< "$DIRname"
	    
	    if [[ ${ENDname[-1]} == *"Cosmol"* ]]; then
		# Then it's a cosmoSLICS run
		# Now get the cosmol ID by splitting on 'Cosmol'
		IFS="Cosmol" read -ra ID <<< "${ENDname[-1]}"
		cosmol=${ID[-1]}  # Only want the last element.
		cosmol_fname=$cosmol
		if [ "$cosmol_fname" != "fid" ]; then # Then it's an integer, make it 2sigfig
		    while [[ ${#cosmol_fname} -lt 2 ]] ; do cosmol_fname="0${cosmol_fname}"; done
		fi
		# Now, which seeds to copy?
		# rsync necessary los
		for l in ${Loop_Array[*]}; do
		    # In case this is a mosaic or noise run, remove extra tags:
		    l=${l%R*} # strip off the "R.."
		    l=${l%n*} # strip off an "n.." 
		    if [ "$l" -lt "26" ]; then
			seed="a"
			los_fname=$l
		    else
			seed="f"
			los_fname=$((l-25))
		    fi
		    # note these lines aren't necessary - access cosmoSLICS cats directly in KiDSMocks_DataGrab.sh
		    #tmp_cosmoSLICS_dataDIR='/data/bengib/Clipping_Pipeline//SLICS_100/'
		    #tmp_cosmoSLICS_DIR=$cosmoSLICS_DIR/${cosmol_fname}_${seed}/GalCat/KV450/GalCat/
		    #if [ ! -d "$tmp_cosmoSLICS_dataDIR" ]; then mkdir -p $tmp_cosmoSLICS_dataDIR; fi
		    #rsync -avz $tmp_cosmoSLICS_DIR/GalCatalog_LOS_cone${los_fname}.fits $tmp_cosmoSLICS_dataDIR/ 
		done
		
	    else
		# Then it's a SLICS run
		# But which SLICS? KiDS-like (KV450 or KiDS1000), or Generic (LSST-like)?
		if [[ "$z" == *"KiDS"* ]]; then
		    # it's KiDS-like SLICS. If z is 'KiDS1000' then it's KiDS1000 mocks. Else it's just plain old KV450.
		    if [ "$z" == "KiDS1000" ]; then
			if [[ "$DIRname" == *"Mosaic"* ]]; then
			    SLICS_KiDS_DIR=${SLICS_Mosaic_DIR}
			else
			    SLICS_KiDS_DIR=${SLICS_KiDS1000_DIR}
			fi
			
			
			# But which redshift bin? JHD's directories for KiDS-1000 are different for each zbin:		       
			source $pipeline_DIR/ShowSumClass/Identify_KiDS1000_zbin.sh $zlo $zhi
			if [ "$zlo" == "0.1" ] && [ "$zhi" == "1.2" ]; then
			    # we need to cycle through the K1000 zbins!
			    bin_names=("bin1" "bin2" "bin3" "bin4" "bin5")
			else
			    # we're accessing a specific bin:
			    bin_names=("$bin_name")
			fi
			
						
		    else
			# Use the KV450 version of SLICS:
			fname_tag="KV450"
			SLICS_KiDS_DIR=${SLICS_KV450_DIR}
			bin_names=("/") # dummy array, just so loop below activates.
		    fi
		    		 
		    ###if [ ! -d "$SLICS_dataDIR" ]; then mkdir -p $SLICS_dataDIR; fi

		    # rsync necessary los
		    if [[ "$DIRname" == *"Mosaic"* ]]; then
			echo "Not copying catalogue - will access directly in KiDSMocks_DataGrab_Mosaic.sh"
		    else
					    
			for b in ${bin_names[*]}; do
			    tmp_SLICS_KiDS_DIR=${SLICS_KiDS_DIR}/$b
			    tmp_SLICS_dataDIR=${SLICS_dataDIR}/$b
			    if [ ! -d "$tmp_SLICS_dataDIR" ]; then mkdir -p $tmp_SLICS_dataDIR; fi
			
			    for l in ${Loop_Array[*]}; do
				if [ "$z" == "KiDS1000" ]; then fname_tag="KiDS1000_${b}"; fi
				rsync -avz $tmp_SLICS_KiDS_DIR/GalCatalog_${fname_tag}_LOS$l.fits $tmp_SLICS_dataDIR/GalCatalog_LOS$l.fits
			    done
			done
		    fi

		    
		elif [[ "$z" == *"LSST"* ]]; then
		    # It's Generic (LSST-like) SLICS
		    if [ ! -d "$SLICS_dataDIR" ]; then mkdir -p $SLICS_dataDIR; fi
		    # rsync necessary los
                    for l in ${Loop_Array[*]}; do rsync -avz $SLICS_Generic_DIR/GalCatalog_LOS$l.fits $SLICS_dataDIR/ ; done
		fi
		
		    
	    fi

	    # The masks
	    if [ ! -d "$G9_Mask_dataDIR" ]; then mkdir -p $G9_Mask_dataDIR; fi
	    if [ ! -d "$W3_Mask_dataDIR" ]; then mkdir -p $W3_Mask_dataDIR; fi
	    rsync -avz $G9_Mask_DIR/G9Mask.*."$sqdeg"deg2.fits $G9_Mask_dataDIR/
	    rsync -avz $W3_Mask_DIR/W3.16bit.*.reg.Now_"$sqdeg"sqdeg.fits $W3_Mask_dataDIR/
		

		
	# Mira Titan
	elif [ "$sqdeg" == "5000" ]; then
		# get the NSIDE from DIRname
		temp=${DIRname#NSIDE}
		NSIDE=${temp%_5000*}
		if [ ! -d "$MiraTitan_dataDIR/NSIDE$NSIDE" ]; then	mkdir -p $MiraTitan_dataDIR/NSIDE$NSIDE; fi 
		# rsync necessary los masks + whole HealPix map
		for l in ${Loop_Array[*]}; do rsync -avz $MiraTitan_DIR/NSIDE$NSIDE/*LOS"$l"_Mask.fits $MiraTitan_dataDIR/NSIDE$NSIDE/ ; done
		rsync -avz $MiraTitan_DIR/NSIDE$NSIDE/*NSIDE$NSIDE*.fits $MiraTitan_dataDIR/NSIDE$NSIDE/
	fi

else	# It's a KiDS Run 

	if [ ! -d "$KiDS_Data_dataDIR" ]; then	mkdir -p $KiDS_Data_dataDIR; fi
	if [ ! -d "$KiDS_Mask_dataDIR" ]; then	mkdir -p $KiDS_Mask_dataDIR; fi
	for f in ${Field[*]}; do
		rsync -avz $KiDS_Data_DIR/KiDS_${f}_reweight_5x5x5_BLIND_inc_m_v2.cat $KiDS_Data_dataDIR/
		rsync -avz $KiDS_Mask_DIR/${f}.16bit.arcmin.AIT.reg2.fits $KiDS_Mask_dataDIR/
	done
fi















