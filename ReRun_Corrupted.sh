#!/bin/bash

# 27/04/2021, B. Giblin
# Launches the clipping pipeline specifically for the cosmologies and redshift bins
# which were listed as corrupted in Joachim's email on 01/04/2021

SS=84.85

# cosmol's 11 and 23 have the same corrupted LOS across zbins3-5:
ZBcut=("0.5-0.7" "0.7-0.9" "0.9-1.2")
SN=(0.273 0.254 0.27)
for j in `seq 0 2`; do
    # ----- cosmol 11 ------
    for los in 16; do
	los2=$((los+25))  # it's seed f (add 25 to LOS)
	sbatch Launch.sh 11 ${SN[$j]} ${ZBcut[$j]} $SS $los2 $los2
    done

    # ----- cosmol 23 ------
    for los in 23; do
	los2=$((los+25))  # it's seed f (add 25 to LOS)
	sbatch Launch.sh 23 ${SN[$j]} ${ZBcut[$j]} $SS $los2 $los2
    done
done

# cosmol 9 has slightly different LOS which are corrupted across zbins3-5
# so have 3 different loops through LOS, depending on the zbin.

# ----- cosmol 9 ------ 
# zbin3
SN="0.273"
ZBcut="0.5-0.7"
for los in 14 18 20 21 23; do
    los2=$((los+25))  # it's seed f (add 25 to LOS)
    sbatch Launch.sh 9 $SN $ZBcut $SS $los2 $los2
done

# zbin4
SN="0.254"
ZBcut="0.7-0.9"
for los in 18 20 21 22; do
    los2=$((los+25))  # it's seed f (add 25 to LOS)                                                                                       
    sbatch Launch.sh 9 $SN $ZBcut $SS $los2 $los2
done

# zbin5
SN="0.27"
ZBcut="0.9-1.2"
for los in 22; do
    los2=$((los+25))  # it's seed f (add 25 to LOS)                                                                                       
    sbatch Launch.sh 9 $SN $ZBcut $SS $los2 $los2
done


