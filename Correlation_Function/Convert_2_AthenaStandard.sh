#!/bin/bash
# 23/09/15, B. M. Giblin, PhD student Edinburgh

# For the KiDS-data, convert coords to RA,DEC (in deg) such that Athena knows the curvature.
# For Sims, multiply X, Y by PS (in deg)
# For both the Sims and the KiDS, remove the <no galaxies> <no columns> from first line
# of unclipped shear file.
# Save in Athena's standard format: X, Y, shear1, shear2, w

# SHOULD BE RUN FROM THE DIRECTORY ABOVE.

pipeline_DIR='/home/bengib/Clipping_SimsLvW/'
data_DIR='/data/bengib/Clipping_SimsLvW/'
source $pipeline_DIR/ShowSumClass/FilterInputArgs.sh $1 $2 $3 $4 $5

indir=$data_DIR/Clipping_K/$DIRname
outdir=$data_DIR/Correlation_Function/"$DIRname"    #_ThBins$ThBins
xy2sky=/usr/bin/xy2sky



if [ "$RUN" == "Sims_Run" ]; then

    # define input/output clipped/unclipped
    combined_name=$name."$gpam"GpAM.LOS$los.SS$SS
    in_rc=$indir/$combined_name.rCLIP_"$sigma"sigma.Xm_Ym_e1c_e2c.asc
    out_rc=$outdir/$combined_name.rCLIP_"$sigma"sigma.ThetaX_ThetaY_e1_e2_w.Std.asc
    orig=$data_DIR/Mass_Recon/$DIRname/$name."$gpam"GpAM.LOS"$los"_Xm_Ym_e1_e2_w.dat
    out_orig=$outdir/$name."$gpam"GpAM.LOS$los.ORIG.ThetaX_ThetaY_e1_e2_w.Std.asc
    
    if [ "$sqdeg" == "60" ] || [ "$sqdeg" == "100" ] || [ "$sqdeg" == "36" ]; then       # You're doing working with a flat sky simulation
		# Need this for low resolution runs.
		echo "MRres is $MRres"
	if [ "$MRres" == "" ] || [ "$MRres" == "-" ] || [ "$MRres" == "5arcs" ]; then
	    angle_conversion=0.001388888888889 #deg/pxl = 5 arcsec / pxl. The default mass recon resolution.
	elif [ "$MRres" == "arcmin" ]; then
	    angle_conversion=0.01666667 # deg/pxl.
	else
	    tmp=$(echo ${MRres%a*}) # strip the number in front of 'arcmin'
	    angle_conversion=$(expr $tmp*0.01666667 | bc)
	fi
	echo "Resolution of mass reconstruction is $angle_conversion deg/pxl"
	awk -v PSm=$angle_conversion '{print $1=$1*PSm,$2=$2*PSm,$3,$4,$5}' < $in_rc > $out_rc
	awk -v PSm=$angle_conversion 'NR>1 {print $1=$1*PSm,$2=$2*PSm,$3,$4,$5}' < $orig > $out_orig

	
    else       # You're working with Mira Titan, and X,Y --> ra,dec must be handled more carefully
        for type in $in_rc $orig
		do
	    if [ "$type" == "$in_rc" ]; then
			Outfile=$out_rc
			awk '{ print $3, $4, $5 }' < $type > $outdir/"$los"_e1_e2_w.asc
	    elif [ "$type" == "$orig" ]; then
			Outfile=$out_orig
			# In the case of the original data file, you need to skip
			# the first line to avid <no. gals> <no. cols>
			awk 'NR>1 { print $3, $4, $5 }' < $type > $outdir/"$los"_e1_e2_w.asc
	    fi

	    # paste the original ra,dec coords onto the e1_e2_w.asc file:
	    paste $data_DIR/Mass_Recon/$DIRname/$name."$gpam"GpAM.LOS"$los"_ra_dec.asc $outdir/"$los"_e1_e2_w.asc > $Outfile

	    # Lastly flip e2, because Kilbinger imagines himself looking down on sky (like God)
	    # ????????? UNSURE; BUT THINK NECESSARY FOR ATHENA, NOT FOR TREECORR
	    awk '{print $1,$2,$3,$4=$4*-1.,$5}' < $Outfile > $outdir/temp$los.asc && mv $outdir/temp$los.asc $Outfile

	    rm -f $outdir/"$los"_e1_e2_w.asc
	done
	
    fi



    

else        # It's a KiDS run
	combined_name=$Field.Blind$Blind.SS$SS

	in_rc=$indir/$combined_name.rCLIP_"$sigma"sigma.Xm_Ym_e1c_e2c.asc 
	orig=$data_DIR/Mass_Recon/$DIRname/$Field.Blind"$Blind"_Xm_Ym_e1_e2_w.dat

	for type in $in_rc $orig 
	do
		if [ "$type" == "$in_rc" ]; then
			Outfile=$outdir/$combined_name.rCLIP_"$sigma"sigma.ThetaX_ThetaY_e1_e2_w.Std.asc
			awk '{ print $1, $2 }' < $type > $outdir/"$Field"_x_y.asc
			awk '{ print $3, $4, $5 }' < $type > $outdir/"$Field"_e1_e2_w.asc

		elif [ "$type" == "$orig" ]; then
			Outfile=$outdir/$Field.Blind$Blind.ORIG.ThetaX_ThetaY_e1_e2_w.Std.asc
			# In the case of the original data file, you need to skip 
			# the first line to avid <no. gals> <no. cols>
			awk 'NR>1 { print $1, $2 }' < $type > $outdir/"$Field"_x_y.asc
			awk 'NR>1 { print $3, $4, $5 }' < $type > $outdir/"$Field"_e1_e2_w.asc
		fi


		# paste the original ra,dec coords onto the e1_e2_w.asc file:
		paste $data_DIR/Mass_Recon/$DIRname/"$Field".Blind"$Blind"_ra_dec.asc $outdir/"$Field"_e1_e2_w.asc > $Outfile		

		# DON'T DO THIS OLD WAY - doing 2 ra --> x,y conversions is fraught with DANGER.
		#$xy2sky -d $data_DIR/Mass_Recon/$DIRname/$Field.16bit.arcmin.AIT.reg2.fits @$outdir/x_y.asc > $outdir/ra_dec_Conversion.asc
		# Take the relevant columns out the conversion file:
		#awk '{print $1, $2}' < $outdir/ra_dec_Conversion.asc > $outdir/ra_dec.asc
		#paste $outdir/ra_dec.asc $outdir/e1_e2_w.asc > $Outfile

		# Lastly flip e2, because Kilbinger imagines himself looking down on sky (like God)
		# When in reality, KiDS data was measured from below.
		awk '{print $1,$2,$3,$4=$4*-1.,$5}' < $Outfile > $outdir/temp$Field.asc && mv $outdir/temp$Field.asc $Outfile

		rm -f $outdir/"$Field"_x_y.asc $outdir/"$Field"_e1_e2_w.asc #$outdir/ra_dec*.asc
	done

fi

