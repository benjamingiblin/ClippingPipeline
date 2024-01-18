# 19/12/2023, B. M. Giblin, Research Fellow, Edn.
# This script shuffles the SLICS KiDS-1000 mosaic LOS
# such that we have N~219 realisations of the K1000 survey
# comprising 18 regions, where no LOS is repeated across
# the 18 regions in any given survey.

import numpy as np
import sys

config = 0
# This is saved as part of the name for the output,
# which is a numpy array shaped ~[Nlos, Nregions],
# which tells you what LOS you should read in for each
# region and survey realisation.
# You could alter the way that the LOS are shuffled
# by altering config in the range [0,216]


# assemble the LOS (74-->292, excluding 198&199)
los = np.arange(74,293)
los = np.delete(los, np.where(los==198))
los = np.delete(los, np.where(los==199))
Nlos = len(los)

Nr = 18
# the number of regions (patches) which make up the K1000 simulated survey

# assemble the original, un-shuffled survey configurations
# i.e. Nlos rows & 18 columns for the regions, where for starters,
# for a given survey congig, every region comes from the same LOS.
orig_survey = np.zeros([ Nlos, Nr ])
for i in range(Nlos):
    orig_survey[i,:] += los[i]

# this creates something structured like this (with 18 cols)
# [ 74, 74, ...., 74 ]
# [ 75, 75, ...., 75 ]
# [ ................ ]
# [ ................ ]
# [ 292, ......, 292 ]

# Now start shuffling:

# if congig=0 then the result will be:
# shift all Region 1 LOS up by 0,
# ..... all Region 2 LOS up by 1, etc.
# and any LOS that go above 292, we shift
# back down to 74, so it's cyclic. 

# BUT... if config=1, then we'll muddle the order of the cols, e.g.,
# in the case of config=1, it will:
# shift all Region 7 LOS up by 0,
# shift all Region 4 LOS up by 1, etc.
# i.e. the Regions which gets shifted up by 1/2/3 etc. are random;
r_order = np.arange(Nr)
if config > 0:
    np.random.seed(config)
    np.random.shuffle( r_order )
    
new_survey = np.copy(orig_survey)
shift = 0
for i in r_order:
    new_survey[:,i] += shift
    shift += 1

# Now fix the elements that get set to 198 or 199
# since these LOS do not exist.
for i in range( Nlos ):
    
    # if 198 exists in a row, shift all elements >=198 by +2
    # this will also fix 199 elements
    if 198 in new_survey[i,:]:
        idx = np.where( new_survey[i,:]>=198 )
        new_survey[i,idx]+=2
    # elif only 199 is in this row, shift by +1.
    elif 199 in new_survey[i,:]:
        idx = np.where( new_survey[i,:]>=199 )
        new_survey[i,idx]+=1
    # elif prevents it from doing BOTH a +1 and +2 shift.

# correct the LOS>292 (the max)
new_survey[ np.where(new_survey>los.max()) ] += ( los.min() - los.max() -1 )

# Now we've got something structured like this!
# (No survey has a repeating LOS across its 18 regions,
#  and each LOS appears only once down the columns, meaning
#  that every LOS and region are used):
# [ 74, 75, 76, ...., 91 ]
# [ 75, 76, 77, ...., 72 ]
# [ ................ ]
# [ ................ ]
# [ 291,292, 74, ......, 89 ]
# [ 292, 74, 75, ......, 90 ]
# (Note that ^this is the result for config 0,
#  for config>0, the ordering is more random!)




# A few checks to make sure the survey configuration passes all requirements:

# 1. No LOS number is repeated for a given regions across all 217 surveys:
for i in range(Nr):
    print("There's %s unqiue LOS for region %s" %( len(np.unique(new_survey[:,i])), i+1 ) )

print( "------------------------------------" )
# 2. Every survey configuration is composed of regions from 18 distinct LOS:
for i in range(Nlos):
    nlos_unique = len( np.unique(new_survey[i,:]) )
    if nlos_unique != Nr:
        print( "Survey realisation at index %s has a LOS repeated!" )

print( "------------------------------------" )

savename = "Shuffled_SLICS-K1000-Mosaic_Configs/Shuffled_Config%s" %config
print( "Saving shuffled survey realisation configuration under this name:" )
print( savename )
np.save(savename, new_survey)
