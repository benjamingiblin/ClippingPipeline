#!/bin/bash
# 05/02/16, B. M. Giblin, PhD student Edinburgh
# Extract data from the KIDS450 FITS catalogue for a given field using ldac

# Use sky2xy to convert the RA,DEC to X,Y in the frame of the mask (which is larger)
# Necessary as Ludo's code uses coords as defined in the MASK frame.

# DO NOT TELL ANYONE WHICH OF THE 'A', 'B' or 'C' you extracted
# ---> Blinding innit. 


#overall_DIR=$PWD
pipeline_DIR='/home/bengib/Clipping_SimsLvW/'
data_DIR='/data/bengib/Clipping_SimsLvW/'
source $pipeline_DIR/ShowSumClass/FilterInputArgs.sh $1 $2 $3 $4 $5


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



DIRECT=$data_DIR/Mass_Recon/$DIRname
ldactoasc=/home/erben/software/theli-1.18.0/bin/Linux_64/ldactoasc_theli
sky2xy=/usr/bin/sky2xy
output=$DIRECT/$Field.Blind"$Blind"_Xm_Ym_e1_e2_w.dat






# Access Catherine's copies of the KiDS-450 data
kids450_datadir='/data/bengib/Clipping_SimsLvW//KiDS450/'
Field_Data=$kids450_datadir/KiDS_"$Field"_reweight_5x5x5_BLIND_inc_m_v2.cat



# Save the RA,DEC in a separate file so we can run sky2xy on it to convert
# these image coords to x,y as defined in the mask frame.
# Make a copy with Z_B in ra,dec file, so we can make the Z_B cut to this catalogue too.
$ldactoasc -i $Field_Data -t OBJECTS -k ALPHA_J2000 DELTA_J2000 -b > $DIRECT/"$Field".Blind"$Blind"_ra_dec.asc
$ldactoasc -i $Field_Data -t OBJECTS -k ALPHA_J2000 DELTA_J2000 Z_B -b > $DIRECT/"$Field".Blind"$Blind"_ra_dec_ZB.asc


# Save z_B in the shear files too - use it to make a cut.
$ldactoasc -i $Field_Data -t OBJECTS -k e1_$Letter e2_$Letter weight_$Letter Z_B -b > $DIRECT/"$Field".Blind"$Blind"_e1_e2_w.asc


$sky2xy -j $PWD/Mass_Recon/KiDS_Fields_Masks/$Field.16bit.arcmin.AIT.reg2.fits @$DIRECT/"$Field".Blind"$Blind"_ra_dec.asc > $DIRECT/"$Field".Blind"$Blind"_x_y_Conversion.asc
# This line saves into the file a bit of nonsense concerning the conversion
# The 5th and 6th columns of this output file are the x and y's we need.
# So reformat:
awk '{print $5, $6}' < $DIRECT/"$Field".Blind"$Blind"_x_y_Conversion.asc > $DIRECT/"$Field".Blind"$Blind"_x_y.asc

paste $DIRECT/"$Field".Blind"$Blind"_x_y.asc $DIRECT/"$Field".Blind"$Blind"_e1_e2_w.asc > $output





######################################### c-corrections #########################################
################################ IMPORTANT 2 MAKE BEFORE Z_b CUTS ###############################

# Calculate the mean of column e1 and columns e2. Apply them as the c-corrections

e1_ccorr=$(echo | awk '{ total += $3 } END { print total/NR }' $output )
echo $e1_ccorr > $DIRECT/"$Field"_e1_ccorr

e2_ccorr=$(echo | awk '{ total += $4 } END { print total/NR }' $output )
echo $e2_ccorr > $DIRECT/"$Field"_e2_ccorr

# Apply the c-corrections
awk -v e1c=$e1_ccorr -v e2c=$e2_ccorr '{print $1, $2, $3-e1c, $4-e2c, $5, $6}' < $output > $DIRECT/temp$Field && mv $DIRECT/temp$Field $output

######################################### c-corrections ######################################### 





if [ $ZBcut != "None" ]; then
		
	echo "Making redshift cut on catalogue. $zlo to $zhi"
	# Save a zB file to investigate n(z)
	awk '{print $6}' < $output > $DIRECT/$Field.Blind"$Blind"_z.dat
	awk -v zl=$zlo -v zh=$zhi '{ if ($6 > zl && $6 <= zh) print $6}' < $output > $DIRECT/$Field.Blind"$Blind"_zCut.dat
	# Make the Z_B cut
	awk -v zl=$zlo -v zh=$zhi '{ if ($6 > zl && $6 <= zh) print $1, $2, $3, $4, $5}' < $output > $DIRECT/temp$Field && mv $DIRECT/temp$Field $output
	awk -v zl=$zlo -v zh=$zhi '{ if ($3 > zl && $3 <= zh) print $1, $2}' < $DIRECT/"$Field".Blind"$Blind"_ra_dec_ZB.asc > $DIRECT/temp$Field && mv $DIRECT/temp$Field $DIRECT/"$Field".Blind"$Blind"_ra_dec.asc
	# ^This ra,dec file will be pasted onto the '_e1_e2_w.asc' file given to athena later. (so it's important)
	# If you're not making a zB cut then no sweat, it's already made at the top.
else
        # Remove the ZB column that was saved
        awk '{print $1, $2, $3, $4, $5}' < $output > $DIRECT/temp$Field && mv $DIRECT/temp$Field $output
fi



# If you want to produce a catalogue with ra,dec instead of X,Y (in mask frame) uncomment the following: 
#awk '{ print $3, $4, $5}' < $output > $DIRECT/temp$Field && mv $DIRECT/temp$Field $DIRECT/"$Field".Blind"$Blind"_e1_e2_w.asc
#paste $DIRECT/"$Field".Blind"$Blind"_ra_dec.asc $DIRECT/"$Field".Blind"$Blind"_e1_e2_w.asc > $DIRECT/"$Field".Blind"$Blind"_ra_dec_e1_e2_w.asc
# The file made here gets used for nothing




# Insert '<no. of lines> 5' at start of catalogue.
# Ludo's code needs no. of galaxies and no. of columns to be listed
# at the start of the input catalogue.
countlines=$(wc -l < $output)
#echo $countlines
echo $countlines' 5' | cat - $output > $DIRECT/temp$Field && mv $DIRECT/temp$Field $output


# Tidy up the Field subdirectory - everythign else is necessary for later.
#rm -f $DIRECT/KiDS_"$Field"_ra_dec_e1_e2_w_Blind"$Blind".asc
rm -f $DIRECT/"$Field".Blind"$Blind"_ra_dec_ZB.asc
rm -f $DIRECT/"$Field".Blind"$Blind"_e1_e2_w.asc
rm -f $DIRECT/"$Field".Blind"$Blind"_x_y.asc
rm -f $DIRECT/"$Field".Blind"$Blind"_x_y_Conversion.asc
# If you want to keep the redshift files, comment the following:
rm -f $DIRECT/$Field.Blind"$Blind"_z.dat
rm -f $DIRECT/$Field.Blind"$Blind"_zCut.dat


