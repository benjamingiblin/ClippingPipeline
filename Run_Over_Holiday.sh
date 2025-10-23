#!/bin/bash

sbatch Launch_JobArray_1a.sh
echo "Just done 1st half of auto-PDFs: I am sleeping for 4 days before killing residual jobs & killing workers..."
sleep 4d
scancel -u bengib
Tree_Correlation_Function/Clean_Workers.sh
sleep 5m

sbatch Launch_JobArray_1b.sh
echo "Just done 2nd half of auto-PDFs: I am sleeping for 4 days before killing residual jobs & killing workers..."
sleep 4d
scancel -u bengib
Tree_Correlation_Function/Clean_Workers.sh
sleep 5m


sbatch Launch_JobArray_2a.sh
echo "Just done 1st half of cross-PDFS: I am sleeping for 6 days before killing residual jobs & killing workers..."
sleep 6d
scancel -u bengib
Tree_Correlation_Function/Clean_Workers.sh
sleep 5m

sbatch Launch_JobArray_2b.sh
echo "Just done 2nd half of cross-PDFS: I am sleeping for 6 days before killing residual jobs & killing workers..."
sleep 6d
scancel -u bengib
Tree_Correlation_Function/Clean_Workers.sh
sleep 5m




