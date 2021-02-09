#!/bin/bash
# Ludo vW's mass mapping code. Edited by Catherine and Benjamin, 24/02/2016

# Required format of input shear catalogue:
# <Number of galaxies> <No. of columns, i.e. 5>
# And then:
# x1, x2, e1, e2, weight
# You also have to be in the Mass_Recon subdir when you run this.

# NOTE TO SELF: This mass map sets the image dimensions to be that of the mask
# Ergo your X and Y coords should be in the mask frame.

module load intel

#overall_DIR=$PWD
pipeline_DIR='/home/bengib/Clipping_SimsLvW/'
data_DIR='/data/bengib/Clipping_SimsLvW/'
source $pipeline_DIR/ShowSumClass/FilterInputArgs.sh $1 $2 $3 $4 $5


echo $Field
# Set the mask and keyword (which distinguishes the temporary output FITS)
if [ "$RUN" == "Sims_Run" ]; then

	if [ "$sqdeg" == "100" ] || [ "$sqdeg" == "60" ] || [ "$sqdeg" == "36" ]; then

		# Find which mask to use:
		getmask=${name##*_} # strip start of name

		if [ "$MRres" == "" ] || [ "$MRres" == "-" ]; then
			MRres="5arcs" # set to default
		fi

		if [ "$getmask" == "Mask" ] || [ "$getmask" == "G9Mask" ] ; then
		    G9mask_datadir='/data/bengib/Clipping_SimsLvW//KiDS450/'
		    mask=$G9mask_datadir/G9Mask.$MRres."$sqdeg"deg2.fits
		else
		    W3mask_datadir='/data/bengib/Clipping_SimsLvW//WMAP_Masks/'
	       	    mask=$W3mask_datadir/W3.16bit.$MRres.reg.Now_"$sqdeg"sqdeg.fits
		fi
	
	else 	# running with Mira Titan
		MT_datadir='/data/bengib/Clipping_SimsLvW//MiraTitan/'
		mask=$MT_datadir/NSIDE"$MRres"/MiraTitan_LOS"$los"_Mask.fits
	fi

	keyword=$los
	combined_name=$data_DIR/Mass_Recon/$DIRname/$name."$gpam"GpAM.LOS"$los"
	echo "The mask is $mask"

else 

	if [ "$MRres" == "1arcmin" ] ; then
		MRres="arcmin" # make sure it matches the name of the mask file
	fi
	mask_startname=$(echo ${Field%_*}) # In case it is a noise run, get the G* part.
	kids450mask_datadir='/data/bengib/Clipping_SimsLvW//KiDS450/'
	mask=$kids450mask_datadir/$mask_startname.16bit.$MRres.AIT.reg2.fits 
	keyword=$Field
	combined_name=$data_DIR/Mass_Recon/$DIRname/$Field.Blind$Blind
fi



# Irrespective of if mask is to be applied, set the end image to the same size as the mask.
# So extract the dimensions of the mask from the header file like so:

stringA=$(fold $mask | grep -a 'NAXIS1')
stringB=$(echo ${stringA##*= }) # Gets rid of the start of the string
								# up to and including '='
nbin1=$(echo ${stringB%/*}) # Gets rid of the end of the string, from '/'

# Do again for the second axis
stringA=$(fold $mask | grep -a 'NAXIS2')
stringB=$(echo ${stringA##*= }) 
nbin2=$(echo ${stringB%/*}) 
#echo $nbin1
#echo $nbin2

inoutdir=$data_DIR/Mass_Recon/$DIRname/		# where the temporary FITS files are saved

for scale in $SS;
do
 
    # Create a symbollic from your mask fits file to newmask.fit.  It's like a copy, but not a physical one...    
	# Does nothing if mask_variable is set to -nomask	
    ln -sf $mask $inoutdir/newmask$keyword.fits

    # The  next code creates galdens.fits, eiso1.fits, eiso2.fits, calib.fits
    # The format of the input ascii file must be:
	# x1, x2, e1, e2, weight

    src/cat2grid_fromasc.a -in \
	"$combined_name"_Xm_Ym_e1_e2_w.dat \
		-directory "$inoutdir" -LOS $keyword -nx $nbin1 -ny $nbin2 \
		-xmax $nbin1 -ymax $nbin2 -xmin 1 -ymin 1

      

	# Mass reconstruction to create the kappa map
    src/massrecon_new.a -directory "$inoutdir" -LOS $keyword -n1 $nbin1 -n2 $nbin2 -gaussian -scale $scale $mask_variable

	
	# Above line makes something called kapparenorm.fits
    # Rename it something useful
    mv $inoutdir/kapparenorm$keyword.fits $combined_name.SS"$scale".Ekappa.fits
    #mv $inoutdir/g1smooth$keyword.fits $combined_name.SS"$scale".g1smooth.fits
    #mv $inoutdir/g2smooth$keyword.fits $combined_name.SS"$scale".g2smooth.fits


    # Now create Bmode mass maps with option "rot" which rotates by 45 degree
    # essentially swap e1->-e2 and e2 ->e1
    # run all the same steps
    
    src/cat2grid_fromasc.a -in \
	"$combined_name"_Xm_Ym_e1_e2_w.dat \
		-directory "$inoutdir" -LOS $keyword -nx $nbin1 -ny $nbin2 \
		-xmax $nbin1 -ymax $nbin2 -xmin 1 -ymin 1 -rot 

    
	src/massrecon_new.a -directory "$inoutdir" -LOS $keyword -n1 $nbin1 -n2 $nbin2 -gaussian -scale $scale $mask_variable


	# Tidy up the produced FITS maps
    mv $inoutdir/kapparenorm$keyword.fits $combined_name.SS"$scale".Bkappa.fits
	

	# Remove all the other FITS files made. Just gets messy.
	rm -f $inoutdir/*$keyword.fits

	
done




