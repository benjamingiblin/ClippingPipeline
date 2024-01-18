#!/bin/bash
# 16/07/2021, B. Giblin, Edinburgh
# Used by the Master_CorrFun_Combine-zbins.sh pipeline.
# This code concatenates shear catalogues together - just two if two paramfiles are entered
# Or, if this is a KiDS1000 run with no tomo, it stitches the 5 cat's from each zbin together.

pipeline_DIR='/home/bengib/Clipping_Pipeline/'
data_DIR='/data/bengib/Clipping_Pipeline/'


# First decide if we are just concatenating shear catalogues,
# or if we are also doing this with the IA catalogues too:
if [[ $IA == *"IA"* ]]; then
    echo "IA detected. Concatenating the IA catalogues as well as the shear."
    file_end_tag=("_Xm_Ym_e1_e2_w.dat" "_IA1_IA2.dat")
else
    file_end_tag=("_Xm_Ym_e1_e2_w.dat")
fi

for fet in ${file_end_tag[*]}; do

    if [[ "$z" == *"KiDS1000"* ]] && [ "$zlo" == "0.1" ] && [ "$zhi" == "1.2" ]; then 
	# Read in the cat from zbin1
	DIRname1="${DIRname/ZBcut0.1-1.2/ZBcut${ZB_array[0]}}"      # alter ZBcut
	DIRname1="${DIRname1/_SN${sigma_e}_/_SN${SN_array[0]}_}"    # alter SN
	name1="${name/SN${sigma_e}_/SN${SN_array[0]}_}"             # alter SN part of filename
	input1=$data_DIR/Mass_Recon/$DIRname1/$name1."$gpam"GpAM.LOS"$los"${fet}
	
	# Now scroll through the remaining bins concatenating
	for zb_idx in `seq 1 4`; do
	    echo "Currently adding in data from zbin $((zb_idx+1))..."

	    tmp_DIRname="${DIRname/ZBcut0.1-1.2/ZBcut${ZB_array[$zb_idx]}}"    # replace ZBcut 
	    tmp_DIRname="${tmp_DIRname/SN${sigma_e}_/SN${SN_array[$zb_idx]}_}" # replace SN
	    tmp_name="${name/SN${sigma_e}_/SN${SN_array[$zb_idx]}_}"           # same for filename
	    tmp_input=$data_DIR/Mass_Recon/$tmp_DIRname/$tmp_name."$gpam"GpAM.LOS"$los"${fet}
	
	    # concatenate:
	    cat $input1 $tmp_input >> $input1
	    rm -f $tmp_input 
	done
	# mv concat file to new directory
	output=$data_DIR/Mass_Recon/$DIRname/$name."$gpam"GpAM.LOS"$los"${fet}
	mv $input1 $output  

    
    else
	# 2 paramfiles have been inputted - concat the two shear cats together.
	# assemble the new (combined-redshift) DIRname & filter inputs:
	source $pipeline_DIR/ShowSumClass/Assemble_Combine-zbin-DIRname.sh $1 $2 $3 $4 $5

	output1=$data_DIR/Mass_Recon/$DIRname1/$name1."$gpam"GpAM.LOS"$los"${fet}
	output2=$data_DIR/Mass_Recon/$DIRname2/$name2."$gpam"GpAM.LOS"$los"${fet}
	
	cat $output1 $output2 > $data_DIR/Mass_Recon/$DIRname/$name1."$gpam"GpAM.LOS"$los"${fet}
	#rm -f $output1 $output2
	# ^this line can cause problems, if a_X_b and a_X_c are running on the same LOS & worker simultaneously
	# so leaving commented out. But this means data products pile up quickly.
	# make sure you don't run too many big jobs in one go, and clean workers often!
    fi

done


