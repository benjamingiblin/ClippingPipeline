#!/bin/bash
# 05/02/16, B. M. Giblin, PhD student Edinburgh
# Extract data from the KIDS450 FITS catalogue for a given field using ldac

# Use sky2xy to convert the RA,DEC to X,Y in the frame of the mask (which is larger)
# Necessary as Ludo's code uses coords as defined in the MASK frame.

# DO NOT TELL ANYONE WHICH OF THE 'A', 'B' or 'C' you extracted
# ---> Blinding innit. 



#source ShowSumClass/FilterInputArgs.sh $1 $2 $3 $4 $5
Field=$1
Blind=$2
SS=$3
mask_name=$4
# Get the end part of the mask name
IFS='.' read -ra ADDR <<< "$mask_name"
file_endname=""
for i in "${ADDR[@]:1}"; do file_endname="$file_endname.$i"; done



######################################### THE KEY 2 BLINDING #########################################
if [ "$Blind" = "2" ]; then
	Letter=A
elif [ "$Blind" = "1" ]; then
	Letter=B
elif [ "$Blind" = "3" ]; then
	Letter=C
else
	echo "You have not entered a valid 1, 2, or 3 for blind."
	exit 1
fi
######################################### THE KEY 2 BLINDING #########################################



DIRECT=$PWD/Mass_Recon/$Field
ldactoasc=/home/cech/software/theli-1.30.0/bin/Linux_64//ldactoasc_theli
sky2xy=
output=$DIRECT/$Field.Blind"$Blind"_Xm_Ym_e1_e2_w"$file_endname".dat






# Access Catherine's copies of the KiDS-450 data
kids450_dir='/home/bengib/KiDS450/'
Field_Data=$kids450_dir/KiDS_"$Field"_reweight_5x5x5_BLIND_PF.cat

# Save the RA,DEC in a separate file so we can run sky2xy on it to convert
# these image coords to x,y as defined in the mask frame.
$ldactoasc -i $Field_Data -t OBJECTS -k ALPHA_J2000 DELTA_J2000 -b > $DIRECT/KiDS_"$Field"_ra_dec_Blind$Blind.asc

# Save z_B now - use it to make a cut.
$ldactoasc -i $Field_Data -t OBJECTS -k e1_$Letter e2_$Letter weight_$Letter Z_B -b > $DIRECT/KiDS_"$Field"_e1_e2_w_Blind$Blind.asc



$sky2xy -j $DIRECT/$mask_name @$DIRECT/KiDS_"$Field"_ra_dec_Blind$Blind.asc > $DIRECT/KiDS_"$Field"_x_y_Blind"$Blind"_Conversion.asc
# This line saves into the file a bit of nonsense concerning the conversion
# The 5th and 6th columns of this output file are the x and y's we need.
# So reformat:
awk '{print $5, $6}' < $DIRECT/KiDS_"$Field"_x_y_Blind"$Blind"_Conversion.asc > $DIRECT/KiDS_"$Field"_x_y_Blind"$Blind".asc

paste $DIRECT/KiDS_"$Field"_x_y_Blind"$Blind".asc $DIRECT/KiDS_"$Field"_e1_e2_w_Blind"$Blind".asc > $output

# Make the Z_B cut
awk '{ if ($6 >= 0.1 && $6 <= 0.9) print $1, $2, $3, $4, $5}' < $output > $DIRECT/temp && mv $DIRECT/temp $output





######################################### c-corrections #########################################

# Calculate the mean of column e1 and columns e2. Apply them as the c-corrections

e1_ccorr=$(echo | awk '{ total += $3 } END { print total/NR }' $output )
echo $e1_ccorr > $DIRECT/"$Field"_e1_ccorr

e2_ccorr=$(echo | awk '{ total += $4 } END { print total/NR }' $output )
echo $e2_ccorr > $DIRECT/"$Field"_e2_ccorr

# Apply the c-corrections
awk -v e1c=$e1_ccorr -v e2c=$e2_ccorr '{print $1, $2, $3-e1c, $4-e2c, $5}' < $output > $DIRECT/temp && mv $DIRECT/temp $output

######################################### c-corrections #########################################




# Insert '<no. of lines> 5' at start of catalogue.
# Ludo's code needs no. of galaxies and no. of columns to be listed
# at the start of the input catalogue.
countlines=$(wc -l < $output)
#echo $countlines
echo $countlines' 5' | cat - $output > $DIRECT/temp && mv $DIRECT/temp $output


# Tidy up the Field subdirectory
rm -f $DIRECT/KiDS_*



