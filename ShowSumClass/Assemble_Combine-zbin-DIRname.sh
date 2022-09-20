#!/bin/bash
# 16/07/2021, B. Giblin, Postdoc, Edinburgh.
# This script is used to assemble a combined-redshift directory name.
# It reads in the DIRname for the first ZBcut, the first ZBcut itself
# and a second ZBcut. It then implants the second ZBcut into the DIRname.
# like this: (*_ZBcutAAA_X_ZBcutBBB_*)

pipeline_DIR='/home/bengib/Clipping_SimsLvW/'
data_DIR='/data/bengib/Clipping_SimsLvW/'

# Do paramfile2 first, since variables like name take their default values
# from paramfile 1 (hence why we're doing it second).
paramfile2=$5
source $pipeline_DIR/ShowSumClass/FilterInputArgs.sh $1 $paramfile2 $3 $4
name2=$name
DIRname2=$DIRname
ZBcut2=$ZBcut
zlo2=$zlo
zhi2=$zhi

# extract the important info form paramfile 2:
paramfile1=$2
source $pipeline_DIR/ShowSumClass/FilterInputArgs.sh $1 $paramfile1 $3 $4
name1=$name
DIRname1=$DIRname
ZBcut1=$ZBcut
zlo1=$zlo
zhi1=$zhi


# assemble the new (combined-redshift) DIRname
delimiter=$ZBcut1
myarray=()
string=$DIRname1$delimiter
while [[ $string ]]; do
    myarray+=( "${string%%"$delimiter"*}" )
    string=${string#*"$delimiter"}
done

# look at the different pieces now saved in myarray
#for value in ${myarray[@]}; do echo "$value "; done
DIRname=${myarray[0]}${ZBcut1}_X_ZBcut${ZBcut2}${myarray[1]}

