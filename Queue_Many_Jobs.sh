#!/bin/bash

# 12/11/2020, B. Giblin
# Launch many clipping/mass mapping jobs, possibly varying cosmology, ZBcut and hence sigma_e.
# Important: Launch.sh must be set up to receive SN & ZBcut params.
# Also make sure the params file in Launch.sh is set to the desired one.

ZBcut=("0.1-0.3" "0.3-0.5" "0.5-0.7" "0.7-0.9" "0.9-1.2")
SN=(0.27 0.258 0.273 0.254 0.27)

# launch jobs for each cosmology and redshift cut:
for i in `seq 0 24`; do
    for j in `seq 0 4`; do
	#echo "Running cosmol $i, ZBcut ${ZBcut[$j]}, sigma_e ${SN[$j]}"
	sbatch Launch.sh $i ${SN[$j]} ${ZBcut[$j]} 1 50
    done
done

# now launch fid
for j in `seq 0 4`; do
    echo "Running cosmol fid, ZBcut ${ZBcut[$j]}, sigma_e ${SN[$j]}" 
    sbatch Launch.sh fid ${SN[$j]} ${ZBcut[$j]} 1 50
done
