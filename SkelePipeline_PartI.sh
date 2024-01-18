#!/bin/bash

# 09/06/16, B. M. Giblin, PhD student Edinburgh

# The skeleton of the clipping pipeline. Executed up to 8 times
# simultaneously by the Master on different Fields/LOS.


# The scripts involved, in order, are:

#  ------------------------- PART I -------------------------
# 0. KiDSMocks_DataGrab.sh or KiDS450_DataGrab.sh (KiDS_Run & 100 sqdeg Sims_Run ONLY) / KiDSMocks_DataGrab_Mosaic.sh
# 1. Sims_DataGrab.py
# 2. Mass_Recon/MassMapLvW.sh 
# 3. get_clipped_shear.py (<Clipping_K>)
#  ------------------------- PART II -------------------------
# 4. Convert_2_AthenaHamana.sh(<Correlation_Function>)
# 5. Athena_CorrFun.sh (<Correlation_Function>)
# 6. plot_CorrFun.py (<Correlation_Function>) <--- Executed from master.

overall_DIR=$PWD
pipeline_DIR='/home/bengib/Clipping_Pipeline/'
data_DIR='/data/bengib/Clipping_Pipeline/'

if [[ "$5" == *"param_files"* ]]; then
    # 2 input paramfiles inputted, this means we're combining shear cats for different zbins
    # assemble the new (combined-redshift) DIRname & filter inputs:
    source $pipeline_DIR/ShowSumClass/Assemble_Combine-zbin-DIRname.sh $1 $2 $3 $4 $5
    paramfile=$paramfile1
else
    # 1 paramfile inputted, just working with shear cat for single zbin
    source $pipeline_DIR/ShowSumClass/FilterInputArgs.sh $1 $2 $3 $4 $5
fi


if [ "$RUN" == "Sims_Run" ]; then
    remove_keyword="*.LOS$los.*"
    remove_keyword2="*LOS"$los"_*"
else
    remove_keyword=$Field"*"
fi

if [[ "$DIRname" == *"Mosaic"* ]]; then
    codeflag="_Mosaic" # affects which code runs: [KiDSMocks_DataGrab(_Mosaic).sh]
fi

# 0. KiDS data and 100 sqdeg mocks require an ldactoasc DataGrab:
if [ "$RUN" == "KiDS_Run" ]; then
	# Is it pure-noise run or actual data run?
	if [[ $Field == *"_NOISE"* ]] && [[ $NOISE == "Y" ]]; then
		python $pipeline_DIR/KiDS450_ShapeNoiseOnly.py $1 $2 $3
	else
		$pipeline_DIR/KiDS450_DataGrab.sh $1 $2 $3
	fi
	
elif [ "$sqdeg" == "100" ]; then

    if [[ "$z" == *"KiDS1000"* ]] && [ "$zlo" == "0.1" ] && [ "$zhi" == "1.2" ]; then
	echo " !!! DOING A KiDS-1000 MOCK RUN w/ NO TOMOGRAPHY - CYCLING THROUGH THE zBINS !!! "
	# Have to run on all tomo bins, then concatenate the resulting cat's.
	paramfile1=$2
	source $pipeline_DIR/ShowSumClass/Identify_KiDS1000_zbin.sh $zlo $zhi
	echo "sigma_e is set to $sigma_e"
	ZB_array=("0.1-0.3" "0.3-0.5" "0.5-0.7" "0.7-0.9" "0.9-1.2")
	SN_array=(0.27 0.258 0.273 0.254 0.27)
	for z_idx in `seq 0 4`; do
	    echo " -x-x-x- On ZBcut${ZB_array[$z_idx]} -x-x-x-x- "
	    tmp_paramfile="${paramfile1/ZBcut0.1-1.2/ZBcut${ZB_array[$z_idx]}}"   # replace the ZBcut label
	    tmp_paramfile="${tmp_paramfile/SN${sigma_e}_/SN${SN_array[$z_idx]}_}"     # replace the SN label	   
	    $pipeline_DIR/KiDSMocks_DataGrab${codeflag}.sh $1 $tmp_paramfile $3 $4
	    python $pipeline_DIR/Sims_DataGrab_MM.py $1 $tmp_paramfile $3 $4  # Add shape noise
	done
	# Now piece the catalogues together:
	source $pipeline_DIR/Mass_Recon/Concatenate_Shear_Cats.sh $1 $2 $3 $4
	python $pipeline_DIR/Mass_Recon/MM_only.py $1 $paramfile $3 $4 
	echo "...Deleting the saved cats in SLICS_100 directory"
	#rm -f $data_DIR/SLICS_100/bin*/GalCatalog_LOS${3}.fits
	# ^ commented out: this could upset other jobs on same worker
	
    else	
	$pipeline_DIR/KiDSMocks_DataGrab${codeflag}.sh $1 $2 $3 $4
	python $pipeline_DIR/Sims_DataGrab_MM.py $1 $2 $3 $4 $5 # Add in the shape noise
    fi
	
    if [[ "$paramfile2" == *"param_files"* ]]; then
	# a 2nd paramfile (different zbin) has been passed in:
	$pipeline_DIR/KiDSMocks_DataGrab${codeflag}.sh $1 $paramfile2 $3 $4
	python $pipeline_DIR/Sims_DataGrab_MM.py $1 $paramfile2 $3 $2 $paramfile1
	
	$pipeline_DIR/Mass_Recon/Concatenate_Shear_Cats.sh $1 $paramfile1 $3 $4 $paramfile2
	echo "Now performing the mass mapping on the concatenated catalogue."
	echo "INPUTS: "
	echo "$1 $2 $3 $4 $5"
	python $pipeline_DIR/Mass_Recon/MM_only.py $1 $2 $3 $4 $5
    fi
    
elif [ "$sqdeg" == "5000" ]; then

	if [ "$ZBcut" == "None" ]; then
		# No ZBcut, work with HealPix maps
		python $pipeline_DIR/MTMocks_DataGrab.py $1 $2 $3 $4 
	else
		# ZBcut means use FITS files, unpack them in Bash, add shape noise in Python.
		$pipeline_DIR/MTMocks_DataGrab.sh $1 $2 $3 $4
		python $pipeline_DIR/MTMocks_DataGrab.py $1 $2 $3 $4
	fi
fi


if [ $? -eq 1 ]; then
	echo "$(date): Data grab failed for params: $1 $2 $3 $4 $5. Exiting SkelePipeline."
	printf "\n$(date): Data grab for $RUN failed for params: $1 $2 $3 $4. Exiting SkelePipeline" > $pipeline_DIR/Error_Reports/Error_Report.txt
   	exit 1
fi

#exit



# Mass mapping and 2D FFT run in to memory issues. Run these calculations for only 
# 1 Field/LOS at a time. Implement a 'lock directory'.
# If LOCKDIR exists, a job is ongoing. The other fields/LOS will pause here until
# LOCKDIR no longer exists.




#LOCKDIR=$data_DIR/LOCKDIR
#locked=1
#while [ $locked = 1 ]
#do
#   sleep 30 # only try to mkdir every 30s. Minimises error message printing.
#   locked=0
#   if mkdir $LOCKDIR ; then
#       echo "Locking succeeded for $1 $2 $3 $4 $5" > $LOCKDIR/Locking_Status
#       echo "Locking succeeded for $1 $2 $3 $4 $5" 		
#   else
#       locked=1
#   fi
#done
	
# 2. MassMapLvW.sh (<Mass_Recon>) 
# You need to change into the Mass_Recon subdir ro run this.
#cd $pipeline_DIR/Mass_Recon
#./MassMapLvW.sh $1 $2 $3 $4 $5

#if [ $? -eq 1 ]; then
#	echo "$(date): MassMapLvW.sh failed for params: $1 $2 $3 $4 $5."
#	printf "\n$(date): MassMapLvW.sh failed for params: $1 $2 $3 $4 $5." > $pipeline_DIR/Error_Reports/Error_Report.txt
#	echo 'Exiting SkelePipeline.'
#	rm -rf $LOCKDIR
#	exit 1
#fi

#cd $pipeline_DIR


#exit
# 3. get_clipped_shear.py (<Clipping_K>) 
#python $pipeline_DIR/Clipping_K/get_clipped_shear.py $1 $2 $3 $4 $5
if [ $? -eq 1 ]; then
	echo "$(date): get_clipped_shear.py failed for params: $1 $2 $3 $4 $5. Exiting SkelePipeline."
	printf "\n$(date): get_clipped_shear.py failed for params: $1 $2 $3 $4 $5. Exiting SkelePipeline." > $pipeline_DIR/Error_Reports/Error_Report.txt
	rm -rf $LOCKDIR
    exit 1
fi




# Unlock
echo "Unlocked!"
rm -rf $LOCKDIR





echo "remove keyword is $remove_keyword for $1 $2 $3 $4"

	
# Remove stuff that Athena doesn't use
removal=$data_DIR/Mass_Recon/$DIRname/"$remove_keyword"kappa.fits
#rm -f $removal


if [[ "$5" == *"param_files"* ]]; then
    # !!! NOTE GIBLIN, if 2 paramfiles inputted, I'm assuming here that you are only calc'ing mass maps."
    # hence, we can delete the concatenated shear cat (otherwise the worker gets full).
    echo "!!! REMOVING SHEAR CAT WITH KEYWORD $remove_keyword2 : comment these lines out if you intend to calc Corr func's !!!"
    removal=$data_DIR/Mass_Recon/$DIRname/"$remove_keyword2".dat
    rm -f $removal
fi






