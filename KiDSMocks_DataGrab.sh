#!/bin/bash
# 26/07/2016 B. M. Giblin, PhD Student, Edinburgh
# This code is for extracting the shear info from the "KiDS-like" mocks
# NB 14/11/2016: Make a ZB cut on the mocks; have to if doing so to the data
# NB will mean that we have < 8.53 gals / arcmin^2 unfortunately. 




#overall_DIR=$PWD
pipeline_DIR='/home/bengib/Clipping_SimsLvW/'
data_DIR='/data/bengib/Clipping_SimsLvW/'
source $pipeline_DIR/ShowSumClass/FilterInputArgs.sh $1 $2 $3 $4 $5

ldactoasc=/home/cech/software/theli-1.30.0/bin/Linux_64//ldactoasc_theli

if [ "$sqdeg" != "100" ]; then
	echo " You are running the data extraction for 100sqdeg mocks, but sqdeg variable is set to $sqdeg"
	echo " Ammend this and try again."
	exit 1
fi

# If you are doing noise cycle (i.e. with cosmoSLICS)
# strip "0n0" to just be "0" and save as new LOS variable name
# to be used ONLY when reading in with ldac
if [[ $DIRname == *"SNCycle"* ]]; then 
	los_fname=${los%n*}
else
	los_fname=${los}
fi


# Get the end part of DIRname
IFS='_' read -ra ENDname <<< "$DIRname"
if [[ ${ENDname[-1]} == *"Cosmol"* ]]; then
    # 'Cosmol' is part of DIRname, it must be a cosmoSLICS run...

    # Now get the cosmol ID by splitting on 'Cosmol'
    IFS="Cosmol" read -ra ID <<< "${ENDname[-1]}"
    cosmol=${ID[-1]}  # Only want the last element.
    cosmol_fname=$cosmol
    if [ "$cosmol_fname" != "fid" ]; then # Then it's an integer, make it 2sigfig
	while [[ ${#cosmol_fname} -lt 2 ]] ; do cosmol_fname="0${cosmol_fname}"; done
    fi

    # If los<26 read in 'a' seeds, otherwise 'f' seeds
    if [ "$los_fname" -lt "26" ]; then
	seed="a"
	los_fname=$los_fname
    else
	seed="f"
	los_fname=$((los_fname-25))
    fi

    # Now use the z keyword to determine if it's KiDS-like or LSST-like cosmoSLICS
    mocks_datadir=/home/jharno/Projects/cosmoSLICS/${cosmol_fname}_${seed}/GalCat/
    if [[ "$z" == *"KiDS"* ]]; then
	# KiDS-like mocks, but is it KV450 or KiDS1000 like?
	if [ "$z" == "KiDS1000" ]; then
	    Tablename="$mocks_datadir/KiDS1000/GalCatalog...."
	    mocks_datadir+="/KiDS1000/GalCat/"
	    keyword_z_phot=z_spec
	    keyword_z_spec=z_spec # redshift cut already made, so only z_spec is saved.

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
	    filename=GalCatalog_KiDS1000_${bin_name}_LOS${los_fname}.fits
	    
	else
	    mocks_datadir+="/KV450/GalCat/"
	    Tablename=GalCatalog
	    # The redshift keywords differ depending on which out of
	    # (SLICS, KiDS-cosmoSLICS, LSST-cosmoSLICS)
	    keyword_z_phot=z_photometric
            keyword_z_spec=z_spectroscopic
	    filename=GalCatalog_LOS_cone${los_fname}.fits
	fi
	
	
    elif [[ "$z" == *"LSST"* ]]; then
	mocks_datadir+="/Generic/"
	Tablename=/disk09/jharno/cosmoSLICS/${cosmol_fname}_${seed}/GalCat/Generic/GalCatalog.dat_LOS${los_fname}
	keyword_z_phot=z_true
	keyword_z_spec=z_true
	filename=GalCatalog_LOS_cone${los_fname}.fits
	
    else
	echo "This script can only accept cosmoSLICS mocks with *zKiDS* or *zLSST* keywords. $z is not allowed"
	exit
    fi
    input=$mocks_datadir/$filename

else
    # 'Cosmol' is not part of name: it's a SLICS run...
    # But which SLICS? KiDS-like (KiDS1000 or KV450) or LSST-like?
    mocks_datadir='/data/bengib/Clipping_SimsLvW//SLICS_100/'
    
    input=$mocks_datadir/GalCatalog_LOS"$los".fits
    # The Table name and redshift keywords differ between KiDS-like (KV450 & KiDS1000) & Generic (LSST-like)
    if [[ "$z" == *"KiDS"* ]]; then
	if [ "$z" == "KiDS1000" ]; then
	    # Infuriatingly for KiDS1000 SLICS, if los>99, Tablename depends on bin (redshift) number
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

	    if [ "$los" -gt 99 ]; then
		Tablename=/disk09/jharno/MockProducts/KV1000/${bin_name}/GalCatalog.dat_LOS${los}
	    else
		Tablename=/disk09/jharno/MockProducts/KV1000/GalCatalog.dat_LOS${los}
	    fi
	    
	    keyword_z_phot=z_true
	    keyword_z_spec=z_true
	else
	    #Annoyingly the Tablename changes for KV450-SLICS depending on LOS:
	    if [ "$los" -gt 215 ]; then
		Tablename=/data/jharno/KiDS/GalCatalog.dat_LOS${los}
	    else
		Tablename=/disk09/jharno/MockProducts/KV450/GalCatalog.dat_LOS${los}
	    fi
	    keyword_z_phot=z_phot
	    keyword_z_spec=z_true
	fi

	
    elif [[ "$z" == *"LSST"* ]]; then
	Tablename=/disk10/jharno/MockProducts/Generic/GalCatalog.dat_LOS${los}
	keyword_z_phot=z_true
	keyword_z_spec=z_true
    fi
    
fi

output=$data_DIR/Mass_Recon/$DIRname/$name."$gpam"GpAM.LOS"$los"_Xm_Ym_e1_e2_w.dat

$ldactoasc -i $input -t $Tablename -k x_arcmin y_arcmin shear1 shear2 $keyword_z_phot $keyword_z_spec > $output


if [ $ZBcut != "None" ]; then
	echo "Making redshift cut on catalogue. $zlo to $zhi"
	#awk '{print $5, $6}' < $output > $data_DIR/Mass_Recon/$DIRname/$name."$gpam"GpAM.LOS"$los"_z.dat
	awk -v zl=$zlo -v zh=$zhi '{ if ($5 >= zl && $5 <= zh) print $5, $6}' < $output > $data_DIR/Mass_Recon/$DIRname/$name."$gpam"GpAM.LOS"$los"_zCut.dat
	# Make the z cut
	awk -v zl=$zlo -v zh=$zhi '{ if ($5 >= zl && $5 <= zh) print $1, $2, $3, $4}' < $output > $data_DIR/Mass_Recon/$DIRname/temp$los_start && mv $data_DIR/Mass_Recon/$DIRname/temp$los_start $output
else
    awk 'NR>6 {print $1, $2, $3, $4}' < $output > $data_DIR/Mass_Recon/$DIRname/temp$los_start && mv $data_DIR/Mass_Recon/$DIRname/temp$los_start $output
fi



# Convert x & y arcmin into pxl values of the simulation.
# Adding shape noise and conversion into MASK frame occurs in Sims_DataGrab.py
PSs=12.908333 # pxls/arcmin in SIMULATION (= 0.00129115558 deg/pxl)
awk -v c=$PSs '{print $1=$1*c, $2=$2*c, $3, $4}' $output > $data_DIR/Mass_Recon/$DIRname/temp$los_start && mv $data_DIR/Mass_Recon/$DIRname/temp$los_start $output






