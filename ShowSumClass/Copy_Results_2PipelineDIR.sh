#!/bin/bash


# 02/11/2017, B. M. Giblin, PhD Student, Edinburgh
# This code is useful if working on a worker of a supercomputer that is not directly connected to your pipeline_DIR.
# This code will only get executed if pipeline_DIR and data_DIR are not the same in Master_CorrFun_ByParts.sh
# This will copy the results produced by the pipeline from data_DIR to pipeline_DIR
# And remove everything left on data_DIR

# End directory
rsync -avz $data_DIR/Tree_Correlation_Function/$DIRname/ThBins$ThBins/NLOS* $pipeline_DIR/Tree_Correlation_Function/$DIRname/ThBins$ThBins/

# The individual CFs produced
rsync -avz $data_DIR/Tree_Correlation_Function/$DIRname/ThBins$ThBins/*.CorrFun.asc $pipeline_DIR/Tree_Correlation_Function/$DIRname/ThBins$ThBins/

# Same again, but for the Correlation_Function directory (used if they used Athena, not TreeCorr)
rsync -avz $data_DIR/Correlation_Function/$DIRname/ThBins$ThBins/NLOS* $pipeline_DIR/Correlation_Function/$DIRname/ThBins$ThBins/
rsync -avz $data_DIR/Correlation_Function/$DIRname/ThBins$ThBins/*.CorrFun.asc $pipeline_DIR/Correlation_Function/$DIRname/ThBins$ThBins/

# Copy Ekappa maps - useful if you want to do cosmology with the 1-pt PDF(K)
rsync -avz $data_DIR/Mass_Recon/$DIRname/*Ekappa.fits $pipeline_DIR/Mass_Recon/$DIRname/
# Copy smoothed shear maps - useful if you want to see what effect of smoothing is on xi_+-
rsync -avz $data_DIR/Mass_Recon/$DIRname/*g*smooth.fits $pipeline_DIR/Mass_Recon/$DIRname/


# Remove everything else. 
#rm -rf $data_DIR/*/$DIRname


if [ "$RUN" == "Sims_Run" ]; then
	# DH10
	if [ "$sqdeg" == "36" ]; then
		echo "Not removing DH10 masks and data from..."
		hostname
		#rm -rf $DH10_dataDIR
		#rm -rf $W3_Mask_dataDIR
		

	# SLICS
	elif [ "$sqdeg" == "100" ]; then
		echo "Not removing SLICS masks and data from..."
		hostname
		#rm -rf $SLICS100_dataDIR
		#rm -rf $SLICS60_dataDIR
		#rm -rf $cosmoSLICS_dataDIR
		#rm -rf $G9_Mask_dataDIR
		#rm -rf $W3_Mask_dataDIR

	# Mira Titan
	elif [ "$sqdeg" == "5000" ]; then	
		echo "Not removing Mira Titan masks and data from..."
		hostname	
		#rm -rf $MiraTitan_dataDIR
	fi

else	# It's a KiDS Run 

	echo "Not removing KiDS masks and data from..."
	hostname
	#rm -rf $KiDS_Mask_dataDIR
	#rm -rf $KiDS_Data_dataDIR

fi









