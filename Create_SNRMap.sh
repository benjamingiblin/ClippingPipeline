#!/bin/bash

overall_DIR=$PWD

# These codes take Field, Blind, SS and N (no. of spin realisations) as params
# as well as the mask name SANS THE directory address of the mask.
Field=$1
echo "$Field"
Blind=$2
SS=$3
N=$4
mask=$5

# Get the end part of the mask name
IFS='.' read -ra ADDR <<< "$mask"
mask_endname=""
for i in "${ADDR[@]:1}"; do mask_endname="$mask_endname.$i"; done


for Field in G12 #G9 G15 #G9 G12 G15 G23 GS
do

	echo $Field
	Timer_File=$overall_DIR/Mass_Recon/$Field/Timer_MakingSNRMap_"$N"Realisations_$mask.txt
	echo "Create SNR Map for $N realisations started at: $(date)" >> $Timer_File
	echo "Create SNR Map for $N realisations started at: $(date)"

	# Redefine the mask to make sure it's right for this field
	mask=$Field$mask_endname

	#./KiDS450_DataGrab_MultiRes.sh $Field $Blind $SS $mask
	echo "Finished producing the KiDS catalogue for this resolution" >> $Timer_File 
	echo "Finished producing the KiDS catalogue for this resolution" 


	counter_max=5 # max no. of processes to run in parallel
	counter=0
	for i in `seq 0 "$N"`;
	do
		#python Mass_Recon/Spin_KiDS_Shear.py $Field $Blind $SS $i $mask &
		counter=$((counter+1));


		# Wait if running 5 processes at once.
		if [ "$counter" == "$counter_max" ]; then 
			echo "Waiting for cores to clear."
			wait 
			counter=0
		fi
	done




	wait


	echo "Finished making $N spin realisations at: $(date)" >> $Timer_File
	echo "Finished making $N spin realisations at: $(date)"






	counter=0
	cd ./Mass_Recon



	./MassMapLvW_SpinKiDS.sh $Field $Blind $SS 0 $mask 
	# Do spin 0 on its own --> this will also make the Ekappa map for the "un-spin" (i.e. actual)
	# shear. The others will then see the Ekappa map exists and won't make it.
	# If you set 5 processes going, each of which tries to make the same Ekappa, it fucks up. 
	for i in `seq 1 "$N"`;
	do
		./MassMapLvW_SpinKiDS.sh $Field $Blind $SS $i $mask &
		counter=$((counter+1));

		# Wait if running 7 processes at once.
		if [ "$counter" == "$counter_max" ]; then 
			echo "Waiting for cores to clear."
			wait 
			counter=0
		fi
	done







	echo "Finished making $N spin mass maps at: $(date)" >> $Timer_File
	echo "Finished making $N spin mass maps at: $(date)"








	wait
	cd $overall_DIR


	python Mass_Recon/Make_SNRMap_From_KappaNSpin.py $Field $Blind $SS $N $mask


	echo "Finished making $N SNR maps at: $(date)" >> $Timer_File
	echo "Finished making $N SNR maps at: $(date)"

	wait

	#rm -f Mass_Recon/$Field/$Field.Blind"$Blind".SS$SS.*kappa.Spin*.fits


done





# Finally make a histogram of the spun KiDS shear: check if the spun-shear histogram looks
# like the B-mode histo for the mocks. If not, it matters that mocks are Gaussian and KiDS are cauchy.
#
#python Mass_Recon/Kappa_Hist.py Sims_Run /disk2/ps1/bengib/Clipping_SimsLvW/param_files/100Sqdeg_SN0.29_Mask_8.53GpAM_X1sigma_SS112_zKiDS 100 219











