#!/bin/bash

if [ $# -ne 8 ]; then
	echo "	Not enough arguments for producing Athena config_tree"
	echo "	ARGUMENTS:"
	echo "	1. KEYWORD TO DISTINGUISH OUTPUT CONFIG_TREE (e.g. Field or los)"
  	echo "	2. GALCAT1"
  	echo "	3. GALCAT2"
  	echo "	4. THMIN"
  	echo "	5. THMAX"
  	echo "	6. NTH (no. of theta bins)"
  	echo "	7. RADEC (0=Cartesian, 1=Spherical)"
  	echo "	8. OATH (0.02 is slow & accurate, 0.05 is quick & dirty)"	
	exit 1
fi

cn=$1
infile1=$2
infile2=$3
thmin=$4
thmax=$5
no_bins=$6
radec=$7
OATH=$8


# Assemble the config file
echo "### -*- sh -*- ###
### Config file for tree code 'athena'
GALCAT1 $infile1 
GALCAT2 $infile2 		# Foreground catalogue
WCORR		1			# 1: shear-shear, 2: shear-position, 4: position-position
SFORMAT		standard			# One of standard, hamana, position
SCOORD_INPUT	deg			# Input catalogue coordinates, {arcsec|arcmin|rad|deg}
SCOORD_OUTPUT	arcmin			# Output coordinates
THMIN           $thmin			# Smallest scale  in units of 'SCOORD_OUTPUT' 
THMAX           $thmax	  		# Largest scale in units of 'SCOORD_OUTPUT'
NTH             $no_bins			# Number of bins
BINTYPE         LOG			# LIN or LOG
RADEC           $radec			# 0: Cartesian, 1: spherical coordinates
OATH		$OATH			# Open angle threshold [rad]
SERROR	        none		# Error type ('none', 'bootstrap', 'jackknife')
NRESAMPLE	50 50			# Number of resampled s" > Correlation_Function/config_files/config_tree$cn
