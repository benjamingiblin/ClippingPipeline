#!/bin/bash
# 02/09/17, B. M. Giblin, PhD Student, Edinburgh
# Extract the data from the Mira Titan mocks with KiDS n(z). Already in FITS file format, no healpix map stuff needed.


pipeline_DIR='/home/bengib/Clipping_SimsLvW/'
data_DIR='/data/bengib/Clipping_SimsLvW/'
source $pipeline_DIR/ShowSumClass/FilterInputArgs.sh $1 $2 $3 $4 $5

ldactoasc=/home/erben/software/theli-1.18.0/bin/Linux_64/ldactoasc_theli
sky2xy=/usr/bin/sky2xy
MT_datadir='/data/bengib/Clipping_SimsLvW//MiraTitan/'

if [ ! -d "$data_DIR/Mass_Recon/$DIRname" ]; then
    mkdir $data_DIR/Mass_Recon/$DIRname
    mkdir $data_DIR/Correlation_Function/$DIRname
fi



# Identify and select which catalogues to use based on the ZBcuts
# The following lists are used to build the files and the tables-within-the-files, to access the data.
tomo_list=("1" "2" "3" "4")
z_list=("0.1" "0.3" "0.5" "0.7" "0.9")
file_list=()

for i in `seq 0 3`; do
	if [ "$zlo" == ${z_list[$i]} ]; then tlo=${tomo_list[$i]}; fi
done

for i in `seq 1 4`; do 
	if [ "$zhi" == ${z_list[$i]} ]; then 
		k=`expr $i - 1` 
		thi=${tomo_list[$k]}
	fi 
done

for i in `seq $tlo $thi`; do
	file_list+=("$MT_datadir/Tomo_Cats/MT_KiDS_tomo$i*")
done
# Have now made a list of the catalogues to include



# Convert the LOS number into ra,dec limits
x=`expr $los % 9`
xp1=`expr $x + 1`
ra_min=$(expr $x*10. | bc)
ra_max=$(expr $xp1*10. | bc)
y=`expr $los / 9`
declinations=("0." "9.59406" "19.47122" "30." "41.81031" "56.44269" "90.")
for i in `seq 0 5`; do
	if [ "$y" -eq "$i" ]; then
		dec_min=${declinations[$i]}
		k=`expr $i + 1`
		dec_max=${declinations[$k]}
	fi
done



DIRECT=$data_DIR/Mass_Recon/$DIRname
output=$DIRECT/$name."$gpam"GpAM.LOS"$los"_Xm_Ym_e1_e2_w.dat
output2=$DIRECT/$name."$gpam"GpAM.LOS"$los"_ra_dec.asc
# Code will break if output files exist already.
rm -f $output
rm -f $output2
#output2=$data_DIR/Correlation_Function/$DIRname/$name."$gpam"GpAM.LOS$los.ORIG.ThetaX_ThetaY_e1_e2_w.Std.asc
for f in ${file_list[*]}; do

	echo "Extracting columns from File $f"
	# Getting the table-name in file
	file_name=$(echo ${f#$MT_datadir/Tomo_Cats/})
	table_name=${file_name::47}"..."

	$ldactoasc -i $f -t MTDir/GalCat/KiDS/$table_name -k ra_arcmin dec_arcmin shear1 shear2 lensfit_w > $output.Tmp
	# Select the right LOS
	awk -v rmin=$ra_min -v rmax=$ra_max -v dmin=$dec_min -v dmax=$dec_max '{ if ($1<=rmax && $1>rmin && $2<=dmax && $2>dmin) print $0}' < $output.Tmp > $DIRECT/$los.Tmp && mv $DIRECT/$los.Tmp $output.Tmp  

	# Now convert ra,dec to X,Y in frame of LOS.
	awk '{print $1, $2}' < $output.Tmp > $DIRECT/$los.ra_dec.Tmp
	awk '{print $3, $4, $5}' < $output.Tmp > $DIRECT/$los.e1_e2_w.Tmp

	$sky2xy -j $MT_datadir/NSIDE$MRres/MiraTitan_LOS${los}_Mask.fits  @$DIRECT/$los.ra_dec.Tmp > $DIRECT/$los.Tmp 
	# The 5th and 6th columns of this file are the x and y's we need.
	awk '{print $5, $6}' < $DIRECT/$los.Tmp > $DIRECT/$los.Tmp2 && mv $DIRECT/$los.Tmp2 $DIRECT/$los.Tmp
	paste $DIRECT/$los.Tmp $DIRECT/$los.e1_e2_w.Tmp > $DIRECT/$los.x_y_e1_e2_w.Tmp

	# If output exists add on the end, otherwise create output
	if [ ! -f $output ]; then
		mv $DIRECT/$los.x_y_e1_e2_w.Tmp $output
		mv $DIRECT/$los.ra_dec.Tmp $output2
	else
		cat $output $DIRECT/$los.x_y_e1_e2_w.Tmp > $DIRECT/$los.Tmp && mv $DIRECT/$los.Tmp $output
		cat $output2 $DIRECT/$los.ra_dec.Tmp > $DIRECT/$los.Tmp && mv $DIRECT/$los.Tmp $output2
	fi
done

# Peg no. gals at start of file.
#countlines=$(wc -l < $output)
#echo $countlines
#echo $countlines' 5' | cat - $output > $DIRECT/$los.Tmp && mv $DIRECT/$los.Tmp $output


# remove unneccesaries.
rm -f $output.Tmp
rm -f $DIRECT/$los.x_y_e1_e2_w.Tmp
rm -f $DIRECT/$los.e1_e2_w.Tmp
rm -f $DIRECT/$los.ra_dec.Tmp
rm -f $DIRECT/$los.Tmp
















