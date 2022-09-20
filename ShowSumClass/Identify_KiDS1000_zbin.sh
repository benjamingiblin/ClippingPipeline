#!/bin/bash

# 18/01/2021, B. Giblin, Postdoc, Edinburgh
# Identify which redshift bin the zlo and zhi input arguments correspond to
# in the case of a KiDS1000-like (cosmo)SLICS run.

zlo=$1
zhi=$2
#echo "zlo and zhi here are $zlo and $zhi"

# Annoying, the redshift cuts have already been made on KiDS1000 mocks, so                                                         
# we're limited to using catatalogues with the KiDS1000 bins:                                                                      
if [ "$zlo" == "0.1" ] && [ "$zhi" == "0.3" ]; then
    bin_name="bin1"
    sigma_e=0.27
    
elif [ "$zlo" == "0.3" ] && [ "$zhi" == "0.5" ]; then
    bin_name="bin2"
    sigma_e=0.258
    
elif [ "$zlo" == "0.5" ] && [ "$zhi" == "0.7" ]; then
    bin_name="bin3"
    sigma_e=0.273
    
elif [ "$zlo" == "0.7" ] && [ "$zhi" == "0.9" ]; then
    bin_name="bin4"
    sigma_e=0.254
    
elif [ "$zlo" == "0.9" ] && [ "$zhi" == "1.2" ]; then
    bin_name="bin5"
    sigma_e=0.27
    
else
    echo "You have set z to KiDS1000. In this case, zlo & zhi must be consecutive numbers from:"
    echo " 0.1, 0.3, 0.5, 0.7, 0.9, 1.2 (i.e. the bins used in KiDS1000 cosmic shear)."
    echo "But you have set zlo and zhi to $zlo, $zhi. Not programmed to deal with this shit. EXITING."
    exit
fi

