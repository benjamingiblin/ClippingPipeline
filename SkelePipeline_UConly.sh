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



source ShowSumClass/FilterInputArgs.sh $1 $2 $3 $4 $5 
overall_DIR=$PWD


if [ "$RUN" == "Sims_Run" ]; then
	remove_keyword="*.LOS$los*"
else
	remove_keyword=$Field"*"
fi




# 0. KiDS data and 100 sqdeg mocks require an ldactoasc DataGrab:
if [ "$RUN" == "KiDS_Run" ]; then
	# Is it pure-noise run or actual data run?
	if [[ $Field == *"_NOISE"* ]] && [[ $NOISE == "Y" ]]; then
		python $overall_DIR/KiDS450_ShapeNoiseOnly.py $1 $2 $3
	else
		$overall_DIR/KiDS450_DataGrab.sh $1 $2 $3
	fi
elif [ "$sqdeg" == "100" ]; then
	$overall_DIR/KiDSMocks_DataGrab.sh $1 $2 $3 $4 
elif [ "$sqdeg" == "5000" ]; then
	python $overall_DIR/MTMocks_DataGrab.py $1 $2 $3 $4 
fi


if [ $? -eq 1 ]; then
	echo "$(date): DataGrab for $RUN failed for params: $1 $2 $3 $4. Exiting SkelePipeline"
	printf "\n$(date): DataGrab for $RUN failed for params: $1 $2 $3 $4. Exiting SkelePipeline" > $overall_DIR/Error_Reports/Error_Report.txt
   	exit 1
fi


	




# 1. Galaxy extraction + masking + adding shape noise
if [ "$RUN" == "Sims_Run" ]; then
	python $overall_DIR/Sims_DataGrab.py $1 $2 $3 $4 $5
fi
	

if [ $? -eq 1 ]; then
	echo "$(date): Sims_DataGrab.py failed for params: $1 $2 $3 $4 $5. Exiting SkelePipeline."
	printf "\n$(date): Sims_DataGrab.py for $RUN failed for params: $1 $2 $3 $4. Exiting SkelePipeline" > $overall_DIR/Error_Reports/Error_Report.txt
   	exit 1
fi






# 2. Convert_2_AthenaHamana_UConly.sh (<Correlation_Function>)
$overall_DIR/Correlation_Function/Convert_2_AthenaStandard_UConly.sh $1 $2 $3 $4 $5
if [ $? -eq 1 ]; then
	echo "$(date): Convert_2_AthenaHamana.sh failed for params: $1 $2 $3 $4 $5. Exiting SkelePipeline."
	printf "\n$(date): Convert_2_AthenaHamana.sh failed for params: $1 $2 $3 $4 $5. Exiting SkelePipeline." > $overall_DIR/Error_Reports/Error_Report.txt
  	exit 1

fi





# Run TreeCorr one at a time, so each run can have the benefit of all cores.

LOCKDIR=$overall_DIR/LOCKDIR
locked=1
while [ $locked = 1 ]
do
   sleep 30 # only try to mkdir every 30s. Minimises error message printing.
   locked=0
   if mkdir $LOCKDIR ; then
       echo "Locking succeeded for $1 $2 $3 $4 $5" > $LOCKDIR/Locking_Status
       echo "Locking succeeded for $1 $2 $3 $4 $5" 		

   else
       locked=1
   fi
done


	

# 4. TreeCorr (<Tree_Correlation_Function>) 
python Tree_Correlation_Function/TreeCorr_CorrFun_UConly.py $1 $2 $3 $4 $5

if [ $? -eq 1 ]; then
	echo "$(date): TreeCorr_CorrFun_UConly.py failed for params: $1 $2 $3 $4 $5."
	printf "\n$(date): TreeCorr_CorrFun_UConly.py failed for params: $1 $2 $3 $4 $5." > $overall_DIR/Error_Reports/Error_Report.txt
	echo 'Exiting SkelePipeline.'
	rm -rf $LOCKDIR
 	exit 1
fi

cd $overall_DIR



# Unlock
echo "Unlocked!"
rm -rf $LOCKDIR





echo "remove keyword is $remove_keyword for $1 $2 $3 $4"

	
# Remove stuff that Athena doesn't use
removal=$overall_DIR/Mass_Recon/$DIRname/"$remove_keyword".dat
rm -f $removal
removal=$overall_DIR/Correlation_Function/$DIRname/"$remove_keyword".asc
rm -f $removal







