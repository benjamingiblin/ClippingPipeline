#!/bin/bash
# 17/11/17, B. M. Giblin, PhD student Edinburgh
# ssh into each worker and remove whatever you left on there

data_DIR='/data/bengib/Clipping_SimsLvW'

echo "There is an exit statement below this line. Remove it to run this code."
echo "!!!! DO NOT DO THIS IF YOU HAVE JOBS RUNNING !!!!"
#exit 0

for i in `seq 1 76`; do
    echo "On worker $i"
    printf -v i "%03d" $i   # make i 3 digit number
    ssh worker$i rm -rf $data_DIR
done
