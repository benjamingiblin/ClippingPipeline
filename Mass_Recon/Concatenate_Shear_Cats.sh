#!/bin/bash
# 16/07/2021, B. Giblin, Edinburgh
# Used by the Master_CorrFun_Combine-zbins.sh pipeline.
# It reads in 2 parameter files, uses them to locate the two shear
# catalogues that have just been made for them, and then concatenates them together.

pipeline_DIR='/home/bengib/Clipping_SimsLvW/'
data_DIR='/data/bengib/Clipping_SimsLvW/'     # If running on a supercomputer, these will be different   

# assemble the new (combined-redshift) DIRname & filter inputs:
source $pipeline_DIR/ShowSumClass/Assemble_Combine-zbin-DIRname.sh $1 $2 $3 $4 $5

output1=$data_DIR/Mass_Recon/$DIRname1/$name1."$gpam"GpAM.LOS"$los"_Xm_Ym_e1_e2_w.dat
output2=$data_DIR/Mass_Recon/$DIRname2/$name2."$gpam"GpAM.LOS"$los"_Xm_Ym_e1_e2_w.dat


# Need to concatenate files but omit the header (first 6 lines) from the second file
output2_tmp=${output2}-tmp
tail -n +7 $output2 > $output2_tmp

cat $output1 $output2_tmp > $data_DIR/Mass_Recon/$DIRname/$name1."$gpam"GpAM.LOS"$los"_Xm_Ym_e1_e2_w.dat
rm -f $output1 $output2 $output2_tmp


if [[ $IA == *"IA"* ]]; then
    echo "IA detected. Concatenating the IA catalogues as well."
    output1_IA=$data_DIR/Mass_Recon/$DIRname1/$name1."$gpam"GpAM.LOS"$los"_IA1_IA2.dat
    output2_IA=$data_DIR/Mass_Recon/$DIRname2/$name2."$gpam"GpAM.LOS"$los"_IA1_IA2.dat

    # miss 1st lines from 2nd file, same as with the shear catalogue
    output2_IA_tmp=${output2_IA}-tmp
    tail -n +7 $output2_IA > $output2_IA_tmp
    
    cat $output1_IA $output2_IA_tmp > $data_DIR/Mass_Recon/$DIRname/$name1."$gpam"GpAM.LOS"$los"_IA1_IA2.dat
    #rm -f $output1_IA $output2_IA $output2_IA_tmp
fi
