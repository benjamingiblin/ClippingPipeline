#!/bin/bash

# 09/06/16, B. M. Giblin, PhD student Edinburgh

# The skeleton of the clipping pipeline. Executed up to 8 times
# simultaneously by the Master on different Fields/LOS.


# The scripts involved, in order, are:

#  ------------------------- PART I -------------------------
# 0. KiDSMocks_DataGrab.sh or KiDS450_DataGrab.sh (KiDS_Run & 100 sqdeg Sims_Run ONLY)
# 1. Sims_DataGrab.py
# 2. runmass.sh (<Mass_Recon>)
# 3. get_clipped_shear.py (<Clipping_K>)
#  ------------------------- PART II -------------------------
# 4. Convert_2_AthenaHamana.sh(<Correlation_Function>)
# 5. Athena_CorrFun.sh (<Correlation_Function>)
# 6. plot_CorrFun.py (<Correlation_Function>) <--- Executed from master.



overall_DIR=$PWD
pipeline_DIR='/home/bengib/Clipping_SimsLvW/'
data_DIR='/data/bengib/Clipping_SimsLvW/'

source $pipeline_DIR/ShowSumClass/FilterInputArgs.sh $1 $2 $3 $4 $5 


##### DECIDE IF YOU WANT TO USE ATHENA OR TREECORR #######
Tree='Y' 				# If ='Y' then it uses TreeCorr, otherwise it uses Athena
##### DECIDE IF YOU WANT TO USE ATHENA OR TREECORR #######


if [ "$RUN" == "Sims_Run" ]; then
	remove_keyword="*.LOS$los.*"
	remove_keyword2="*LOS"$los"_*"
else
	remove_keyword=$Field"*"
fi






# 4. Convert_2_AthenaHamana.sh (<Correlation_Function>)
$pipeline_DIR/Correlation_Function/Convert_2_AthenaStandard.sh $1 $2 $3 $4 $5
if [ $? -eq 1 ]; then
	echo "$(date): Convert_2_AthenaHamana.sh failed for params: $1 $2 $3 $4 $5. Exiting SkelePipeline."
	printf "\n$(date): Convert_2_AthenaHamana.sh failed for params: $1 $2 $3 $4 $5. Exiting SkelePipeline." > $pipeline_DIR/Error_Reports/Error_Report.txt
  	exit 1

fi





if [ "$Tree" == "Y" ]; then
#        module load anaconda/2     # required to import yaml 
	# Only do TreeCorr one at a time as it is parallelised... require a locking directory

	LOCKDIR=$data_DIR/LOCKDIR2
	locked=1
	while [ $locked = 1 ]
	do
  		sleep 30 # only try to mkdir every 30s. Minimises error message printing.
   		locked=0
   		if mkdir $LOCKDIR ; then
       		echo "Locking TreeCorr succeeded for $1 $2 $3 $4 $5" > $LOCKDIR/Locking_Status
       		echo "Locking TreeCorr succeeded for $1 $2 $3 $4 $5" 		
   		else
       		locked=1
   		fi
	done

	# 5. TreeCorr_CorrFun.py (<Tree_Correlation_Function>)
	python $pipeline_DIR/Tree_Correlation_Function/TreeCorr_CorrFun.py $1 $2 $3 $4 $5
	if [ $? -eq 1 ]; then
	 	echo "$(date): TreeCorr_CorrFun.py failed for params: $1 $2 $3 $4 $5. Exiting SkelePipeline."
		printf "\n$(date): TreeCorr_CorrFun.py failed for params: $1 $2 $3 $4 $5. Exiting SkelePipeline." > $pipeline_DIR/Error_Reports/Error_Report.txt
		# Unlock
		rm -rf $LOCKDIR
	  	exit 1
	fi
	
	# Unlock
	echo "Unlocked!"
	rm -rf $LOCKDIR


else

	# 5. Athena_CorrFun.sh (<Correlation_Function>)
	$pipeline_DIR/Correlation_Function/Athena_CorrFun.sh $1 $2 $3 $4 $5
	if [ $? -eq 1 ]; then
	 	echo "$(date): Athena_CorrFun.sh failed for params: $1 $2 $3 $4 $5. Exiting SkelePipeline."
		printf "\n$(date): Athena_CorrFun.sh failed for params: $1 $2 $3 $4 $5. Exiting SkelePipeline." > $pipeline_DIR/Error_Reports/Error_Report.txt
	  	exit 1
	fi
fi





echo "remove keyword is $remove_keyword for $1 $2 $3 $4"
echo "remove keyword2 is $remove_keyword for $1 $2 $3 $4"
# Remove stuff that Athena doesn't use

removal=$data_DIR/Mass_Recon/$DIRname/"$remove_keyword"kappa.fits
#rm -f $removal
removal=$data_DIR/Mass_Recon/$DIRname/"$remove_keyword2".dat
rm -f $removal



removal2=$data_DIR/Clipping_K/$DIRname/$remove_keyword.fits
#rm -f $removal2
removal2b=$data_DIR/Clipping_K/$DIRname/$remove_keyword.asc
rm -f $removal2b
removal2c=$data_DIR/Correlation_Function/$DIRname/$remove_keyword.Std.asc
#rm -f $removal2c




