#!/bin/bash

# 01/06/2017, B. M. Giblin, PhD Student Edinburgh
# The following code is used to alter the following within the relevant codes of the pipeline:
# 		1. The directory for KiDS data
#		2. The directory for KiDS masks
#		3. The directory for SLICS 100 deg^2 mocks (KiDS-like & Generic (LSST-like))
#		4. The directories for the W3/G9 masks to apply to SLICS 
#		5. The directory for Mira Titan (REDUNDANT)
#		6. The directory for the DH10 mocks (REDUNDANT)

#		7. The executable sky2xy & xy2sky
#		8. The executable ldactoasc
#		9. The executable for Athena


echo " An exit statement is in place just after this line, to prevent Install_Pipeline.sh being ran accidentally."
echo " Comment this out to continue"
exit 0


# !!! ATTENTION !!!
# Before running this code to install the pipeline, you need to do the following:
# 1. move the Install_Pipeline directory upwards one level
# 2. change the following directories to reflect where your KiDS-450 data, SLICS and Mira Titan mocks are,
#    if you're running on cuillin, leave these alone.

# If you are running this on cuillin, you'll want to set:
# data_DIR=/data/<username>/Clipping_Pipeline/
# and set pipeline_DIR to the directory you have cloned the repo into.

username=$(whoami)
data_DIR=/data/$username/Clipping_Pipeline/
pipeline_DIR=/home/$username/Clipping_Pipeline/

# Define where the data/masks/mocks are if you are saving in pipeline_DIR 
KiDS_Data_DIR=/home/$username/KiDS450/
KiDS_Mask_DIR=/home/$username/KiDS450/
	
SLICS_KV450_DIR=/disk10/jharno/MockProducts/KV450/
SLICS_KiDS1000_DIR=/disk10/jharno/MockProducts/KiDS1000/
SLICS_Generic_DIR=/disk10/jharno/MockProducts/Generic/

G9_Mask_DIR=/home/$username/KiDS450/
W3_Mask_DIR=/home/bengib/WMAP_Masks/ # Leave this pointing to Ben's directory
	
MiraTitan_DIR=/home/$username/MiraTitan/

DH10_DIR=/home/$username/DH10_Mocks/FaLCoNS/

# Define where data/masks/mocks are if you are saving in data_DIR 
KiDS_Data_dataDIR=$data_DIR/KiDS450/
KiDS_Mask_dataDIR=$data_DIR/KiDS450/

SLICS_dataDIR=$data_DIR/SLICS_100/

G9_Mask_dataDIR=$data_DIR/KiDS450/
W3_Mask_dataDIR=$data_DIR/WMAP_Masks/

MiraTitan_dataDIR=$data_DIR/MiraTitan/

DH10_dataDIR=$data_DIR/DH10_Mocks/FaLCoNS/



# !!! ATTENTION !!!
# Executables: sky2xy, xy2sky, ldactoasc and athena NEED TO BE IN $PATH for the following to work.
# If they're not, please add them to $PATH via your bashrc
sky2xy="$(which sky2xy)"
xy2sky="$(which xy2sky)"
ldactoasc=/home/cech/software/theli-1.30.0/bin/Linux_64//ldactoasc_theli #"$(which ldactoasc)"
athena="$(which athena)"


# Begin editing codes.

# 0. pipeline_ and data_DIR
find $pipeline_DIR -type f -exec sed -i "s#data_DIR=.*#data_DIR=\'$data_DIR\'#" {} + 
find $pipeline_DIR -type f -exec sed -i "s#pipeline_DIR=.*#pipeline_DIR=\'$pipeline_DIR\'#" {} +

# 1. IF USING pipeline_DIR
# KiDS-450 data and mask
find $pipeline_DIR -type f -exec sed -i "s#kids450_dir=.*#kids450_dir=\'$KiDS_Data_DIR\'#" {} +
find $pipeline_DIR -type f -exec sed -i "s#kids450maskdir=.*#kids450maskdir=\'$KiDS_Mask_DIR\'#" {} +
# SLICS
find $pipeline_DIR -type f -exec sed -i "s#SLICS_KV450_DIR=.*#SLICS_KV450_DIR=\'$SLICS_KV450_DIR\'#" {} +
find $pipeline_DIR -type f -exec sed -i "s#SLICS_KiDS1000_DIR=.*#SLICS_KiDS1000_DIR=\'$SLICS_KiDS1000_DIR\'#" {} +
find $pipeline_DIR -type f -exec sed -i "s#SLICS_Generic_DIR=.*#SLICS_Generic_DIR=\'$SLICS_Generic_DIR\'#" {} +

# Masks to apply to SLICS
find $pipeline_DIR -type f -exec sed -i "s#G9maskdir=.*#G9maskdir=\'$G9_Mask_DIR\'#" {} +
find $pipeline_DIR -type f -exec sed -i "s#W3maskdir=.*#W3maskdir=\'$W3_Mask_DIR\'#" {} +
# Mira Titan 
find $pipeline_DIR -type f -exec sed -i "s#MTdir=.*#MTdir=\'$MiraTitan_DIR\'#" {} +
# DH10
find $pipeline_DIR -type f -exec sed -i "s#DH10_dir=.*#DH10_dir=\'$DH10_DIR\'#" {} +


# 2. IF USING data_DIR

# KiDS-data
find $pipeline_DIR -type f -exec sed -i "s#kids450_datadir=.*#kids450_datadir=\'$KiDS_Data_dataDIR\'#" {} +
find $pipeline_DIR -type f -exec sed -i "s#kids450mask_datadir=.*#kids450mask_datadir=\'$KiDS_Mask_dataDIR\'#" {} +

# SLICS
find $pipeline_DIR -type f -exec sed -i "s#SLICS_dataDIR=.*#SLICS_dataDIR=\'$SLICS_dataDIR\'#" {} +

# Masks to apply to SLICS
find $pipeline_DIR -type f -exec sed -i "s#G9mask_datadir=.*#G9mask_datadir=\'$G9_Mask_dataDIR\'#" {} +
find $pipeline_DIR -type f -exec sed -i "s#W3mask_datadir=.*#W3mask_datadir=\'$W3_Mask_dataDIR\'#" {} +
# Mira Titan 
find $pipeline_DIR -type f -exec sed -i "s#MT_datadir=.*#MT_datadir=\'$MiraTitan_dataDIR\'#" {} +
# DH10
find $pipeline_DIR -type f -exec sed -i "s#DH10_datadir=.*#DH10_datadir=\'$DH10_dataDIR\'#" {} +

# Executables
find $pipeline_DIR -type f -exec sed -i "s#ldactoasc=.*#ldactoasc=$ldactoasc#" {} +
find $pipeline_DIR -type f -exec sed -i "s#sky2xy=.*#sky2xy=$sky2xy#" {} +
find $pipeline_DIR -type f -exec sed -i "s#xy2sky=.*#xy2sky=$xy2sky#" {} +
find $pipeline_DIR -type f -exec sed -i "s#athena=.*#athena=$athena#" {} +




# Some things to do with the mass reconstruction fortran codes
# change line that directs to the cookbook.f90
cd $pipeline_DIR/Mass_Recon/src/
# no longer necessary
#address=$(find `pwd` -name "cookbook.f90" )
#find . -type f -exec sed -i "s#/home/cech/cookbook.f90#$address#" {} +

# compile mass recon codes
module load compiler mkl
#####ifort cat2grid_fromasc.f90 -o cat2grid_fromasc.a -L/usr/local/share/cfitsio/ -lcfitsio
ifort cat2grid_fromasc.f90 -o cat2grid_fromasc.a -L/usr/local/share/cfitsio/ -lcfitsio -O0 -CB 

ifort massrecon_new.f90 -o massrecon_new.a -L/usr/local/share/cfitsio/ -lcfitsio -L/roe/intel/Compiler/11.1/072/mkl/lib/em64t/ -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -qopenmp -lpthread -L/usr/lib -lm -lfftw3

cd ../../../Install_Pipeline/
