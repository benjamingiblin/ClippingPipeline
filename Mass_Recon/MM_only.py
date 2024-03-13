# 26/07/2023 B. M. Giblin, Research Fellow, Edinburgh
# Grab the concatenated ellipticity catalogues and perform mass mapping

# NEW UPDATE 08/06/2023: Perform mass mapping with lenspack...:
from lenspack.utils import bin2d
from lenspack.image.inversion import ks93
from scipy.ndimage import gaussian_filter

import os.path, sys
import numpy as np
from astropy.io import fits
import random
from subprocess import call
from os import getcwd
import math
import glob

pipeline_DIR='/home/bengib/Clipping_Pipeline/'
data_DIR='/data/bengib/Clipping_Pipeline/'

classdir = pipeline_DIR + "/ShowSumClass"
sys.path.insert(0, classdir) # add directory in which classes & functions 
							 # are defined to the python path
from ClassWarfare import Filter_Input, Format_2Darray
from FunkShins import interpolate2D, Mask_Shear, Lower_Res_Mask, Combine_zbin_DIRname

if len(sys.argv)-1 > 4 and "param_files" in sys.argv[-1]:
	# 2  paramfiles have given, it's a combination/cross redshift problem (change input/output DIRname)
	name, gpam, DIRname, SS, sigma, SN, mask, z, PS, sqdeg, zlo, zhi, ThBins, OATH, los, los_end = Combine_zbin_DIRname( sys.argv )
else:
	# Note this shouldn't be activ'd, as this code only runs for X-redshift runs.
	variable = Filter_Input(sys.argv)
	variable.Filter()
	name, gpam, DIRname, SS, sigma, SN, mask, z, PS, sqdeg, zlo, zhi, ThBins, OATH, los, los_end = variable.Unpack_Sims()

RUN = sys.argv[1]
#print(name, gpam, DIRname, SS, sigma, SN, mask, z, PS, sqdeg, los, los_end)

print(" READING IN A CONCATENATED ELLIPTICITY CATALOGUE TO PERFORM THE MASS MAPPING ")
#print("Reading in the %s sqdeg shear catalogue for LOS %s" %(sqdeg, los))
# Pxls already in mask frame


# Get the size of the map
Prepend = ''.join([ i for i in list(DIRname)[0:5] if i in ['M','R','r','e','s'] ])
if Prepend == 'MRres':
	MRres = DIRname.split('_')[0].split('MRres')[1]
	if 'arcmin' in MRres:
		result = MRres.split('arcmin')[0]
		PSm = float(result)/60. # Pxl scale of the mask, deg/pxl
	else:
		result = MRres.split('arcs')[0]
		PSm = float(result)/3600.
else:
	MRres= '5arcs'
	PSm = PSm_5arcs

new_sizeX = int( math.ceil( (np.sqrt(float(sqdeg)) / PSm) /100. ) )*100
new_sizeY = new_sizeX


if 'Mosaic' in DIRname:
	# Read in mask and use this as the map dimensions
	Rfactor = int( los.split('R')[-1].split('n')[0] )
	from astropy.wcs import WCS
	maskdir = '/home/bengib/KiDS1000_Data/Masks_SLICS_Regions'
	f = fits.open('%s/Mask_KiDS1000_R%s_Filter12.5_%s.fits' %(maskdir,Rfactor,MRres)) 
	# doesn't matter which R mask you use actually, all  are rotated to equator.
	w = WCS( f[0].header )
	new_sizeX = f[0].data.shape[1]
	new_sizeY = f[0].data.shape[0]
	mask_pxls = np.nonzero( f[0].data )

	# Thisis the simple flat-sky approach, agrees w/ above method at 1-2% level
	# NB: X,Y are in deg, as is PSm, so no /60. factor.        
	#X = X/PSm
	#Y = Y/PSm
        

X,Y,e1,e2,Weight = np.loadtxt('%s/Mass_Recon/%s/%s.%sGpAM.LOS%s_Xm_Ym_e1_e2_w.dat'%(data_DIR, DIRname, name, gpam, los),
                              usecols=(0,1,2,3,4), unpack=True)


#sys.exit()
# ----- PERFORM THE MASS MAPPING -----
print(" Performing the mass mapping with lenspack, smoothing scale %s pxls " %SS )
# project ellip. onto grid:
if 'Mosaic' in DIRname:
	# not sure the bin2d stat works well with masks! manually 2d yourself!
	from FunkShins import MeanQ_VS_XY
	e1map,_,_,bad_pxls = MeanQ_VS_XY(e1, Weight, np.ones_like(e1), X,Y,(new_sizeY,new_sizeX))
	e2map,_,_,bad_pxls = MeanQ_VS_XY(e2, Weight, np.ones_like(e1), X,Y,(new_sizeY,new_sizeX))
else:
	e1map, e2map = bin2d(X, Y, v=(e1,e2), w=Weight, npix=(new_sizeX,new_sizeY), extent=None)
	# use extent to set mask boundaries (if mosaic mocks) 

e1map_sm = gaussian_filter(e1map, sigma=float(SS))        # smooth
e2map_sm = gaussian_filter(e2map, sigma=float(SS))        # smooth

# PAD up by factor of 2:
pad_size = 2 * max(new_sizeX,new_sizeY)

# Find smallest power of 2 that gives size > map size
#sm = np.where(2**np.arange(16) > max(new_sizeX,new_sizeY) )[0][0]
#pad_size = 2**sm

pad_X = int((pad_size - new_sizeX)/2) # num. of rows to add on each side
pad_Y = int((pad_size - new_sizeY)/2) # num. cols to add on each side
pad_e1map_sm = np.pad( e1map_sm, ((pad_Y, pad_Y),(pad_X, pad_X)), mode='constant')  # pad e1 map with zeros   
pad_e2map_sm = np.pad( e2map_sm, ((pad_Y, pad_Y),(pad_X, pad_X)), mode='constant')

# Mass mapping:
pad_kappaE,pad_kappaB = ks93(pad_e1map_sm, pad_e2map_sm)                 # padded mass map
kappaE = pad_kappaE[ pad_Y:(pad_Y+new_sizeY), pad_X:(pad_X+new_sizeX) ]  # remove padding
kappaB = pad_kappaB[ pad_Y:(pad_Y+new_sizeY), pad_X:(pad_X+new_sizeX) ]
if 'Mosaic' in DIRname:
	# kill the masked_pxls
	kappaE[mask_pxls[0], mask_pxls[1]] = 0.
	kappaB[mask_pxls[0], mask_pxls[1]] = 0.

# save:
savename_map = '%s/Mass_Recon/%s/%s.%sGpAM.LOS%s.SS%s' %(data_DIR,DIRname, name,gpam,los,SS)
np.save( savename_map+'.Ekappa', kappaE )
np.save( savename_map+'.Bkappa', kappaB )







