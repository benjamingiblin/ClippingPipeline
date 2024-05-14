#!/bin/bash
# 18/03/2024
# Check if data disk exists on a node.

workername=$1
data_DIR='/data/'
home_DIR='/home/bengib/Clipping_Pipeline/'
if [ -d $data_DIR ]; then echo $workername >> ${home_DIR}/nodelist.txt; fi
