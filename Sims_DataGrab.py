# 09/09/16 B. M. Giblin, PhD student, Edinburgh
# Grab the shear data for the simulated LOS
# Mask it (reformatted W3 mask) and add shape noise if specified.
# Can deal with both the 60 sqdeg ("Sims_Run") and 100 sqdeg KiDS-like ("KiDSSims_Run") mocks


import os.path, sys
import numpy as np
from astropy.io import fits
import random
from subprocess import call
from os import getcwd
import math
import glob

pipeline_DIR='/home/bengib/Clipping_SimsLvW/'
data_DIR='/data/bengib/Clipping_SimsLvW/'

classdir = pipeline_DIR + "/ShowSumClass"
sys.path.insert(0, classdir) # add directory in which classes & functions 
							 # are defined to the python path
from ClassWarfare import Filter_Input, Format_2Darray
from FunkShins import interpolate2D, Mask_Shear, Lower_Res_Mask


variable = Filter_Input(sys.argv)
variable.Filter()
name, gpam, DIRname, SS, sigma, SN, mask, z, PS, sqdeg, zlo, zhi, ThBins, OATH, los, los_end = variable.Unpack_Sims()
RUN = sys.argv[1]

#print(name, gpam, DIRname, SS, sigma, SN, mask, z, PS, sqdeg, los, los_end)

PSm_5arcs = 1.388888888889e-03 # Pxl scale of the full res mask, deg/pxl
# Determine the resolution of the mask you will use in mass reconstruction
# First get the MRres
Prepend = ''.join([ i for i in list(DIRname)[0:5] if i in ['M','R','r','e','s'] ])
if Prepend == 'MRres':
        MRres = DIRname.split('_')[0].split('MRres')[1]
	#result = ''.join([i for i in MRres if i.isdigit()])
        
        if 'arcmin' in MRres:  # pxl scale expressed in arcmins
                result = MRres.split('arcmin')[0]
                PSm = float(result)/60. # Pxl scale of the mask, deg/pxl
                
        else:                  # pxl scale expressed in arcsec
                result = MRres.split('arcs')[0]
                PSm = float(result)/3600.
                
else:
        MRres= '5arcs'
        PSm = PSm_5arcs 





# If the mask does not currently exist, it needs to be made. This only works for the W3 mask currently.
# In order to be made it needs to know how long to make the sides: = (len of mock in deg / PSm) rounded to nearest 100
new_sizeX = int( math.ceil( (np.sqrt(float(sqdeg)) / PSm) /100. ) )*100
new_sizeY = new_sizeX

# Also needs to know how big the mask should be if at FULL res of 5arcs
Full_sizeX = int( math.ceil( (np.sqrt(float(sqdeg)) / PSm_5arcs) /100. ) )*100
Full_sizeY = Full_sizeX


if name.split('_')[1] == 'Mask' or name.split('_')[1] == 'G9Mask':
	G9mask_datadir='/data/bengib/Clipping_SimsLvW//KiDS450/'
	mask_filename = '%s/G9Mask.%s.%sdeg2.fits'%(G9mask_datadir, MRres, sqdeg)
elif name.split('_')[1] == 'W3Mask':
	W3mask_datadir='/data/bengib/Clipping_SimsLvW//WMAP_Masks/'
	mask_filename = '%s/W3.16bit.%s.reg.Now_%ssqdeg.fits'%(W3mask_datadir, MRres, sqdeg)

	if os.path.exists('%s/W3.16bit.%s.reg.Now_%ssqdeg.fits'%(W3mask_datadir, MRres, sqdeg)) != True :
		Mask = fits.open('%s/W3.16bit.5arcs.reg.fits'%W3mask_datadir)
		Mask[0].data = Format_2Darray(Mask[0].data).Manipulate_W3(Full_sizeX, Full_sizeY, 1.)

		if MRres != '5arcs':
			# If resolution is < 5arcs, reduce the res of the mask just produced:
			print("This was activated. Making new mask")
			Mask[0].data, masked = Lower_Res_Mask( Mask[0].data, new_sizeX, new_sizeY)
		Mask.writeto('%s/W3.16bit.%s.reg.Now_%ssqdeg.fits'%(W3mask_datadir, MRres, sqdeg), output_verify='ignore', clobber=True)


def Generate_Unitary_Shape_Noise(seedy, SN_level):
	# Generate Gaussian shape noise bounded by [-1,+1]
	np.random.seed(seedy)
	e_rng = np.random.normal(0., SN_level, len(e1_temp))
	idx = np.where(abs(e_rng)>1.)[0]						# Find the trouble makers >1, <-1
	for i in idx:											# Fix the trouble-makers
		newseed=i
		np.random.seed(newseed)
		new_e = np.random.normal(0., SN_level)
		while abs(new_e) > 1:
			newseed +=1
			np.random.seed(newseed)
			new_e = np.random.normal(0., SN_level)
		e_rng[i] = new_e
	return e_rng


if sqdeg==60:

	mocks60_datadir='/data/bengib/Clipping_SimsLvW//SLICS_60/'
	mocks60_datadir += 'z_%s'%z
	hdulist1 = fits.open('%s/%sgamma1_weight.dat_LOS%s.fits'%(mocks60_datadir,z,los))
	hdulist2 = fits.open('%s/%sgamma2_weight.dat_LOS%s.fits'%(mocks60_datadir,z,los))
	grid1 = hdulist1[0].data
	grid2 = hdulist2[0].data


	# Convert gal/sq arcmin into a number of pxls (galaxies) to extract
	size = int(float(gpam)*216000.) # Multiply by no. of square arcmin in image
	# The number you end up having in your cat will be less if you do masking


	xpxl_len = len(grid1[0,:])
	ypxl_len = len(grid1[:,0]) # sims are square, so these are the same.

	# Randomly pick pxls and interpolate to get shear:
	# Randomly pick THE SAME PXLS every time for a given LOS,
	# BUT different pxls for each successive LOS.
	np.random.seed(int(los))
	X = np.random.uniform(0,xpxl_len-1, size) 
	np.random.seed(2*int(los))
	Y = np.random.uniform(0,ypxl_len-1, size) 

	e1_temp = interpolate2D(X, Y, grid1) 
	e2_temp = interpolate2D(X, Y, grid2)

elif sqdeg==36:
	num_rt_cats = 5		# the number of ray tracing catalogues per realisation
	DH10_datadir='/data/bengib/Clipping_SimsLvW//DH10_Mocks/FaLCoNS/'
	# get the cosmol number
	cosmol = int(DIRname.split('_Cosmol')[-1])

	# get the los number + noise cycle number if applicable
	intlos = int(los.split('n')[0])		
	# get the realisation number
	# The following 3 lines I coded in, and then decided must be wrong. Incorrectly gets realisation number for multiples of 5.
	
	#if intlos > 0 and intlos % num_rt_cats == 0:
	#	realisation_number = intlos/num_rt_cats
	#else:
	realisation_number = (intlos/num_rt_cats) +1		# realisation numbers index from 1.
	
	# get the catalogue number
	catalogue_number = intlos % num_rt_cats	
	allcatalogues = glob.glob('%s/Cosmol%s/N*_L*_Om*_s8*_w-1.00-%03d/rt/gal/*/FaLCoNS_with_KiDS_nofz_0.5_z_0.9.cat'%(DH10_datadir,cosmol,realisation_number))
	catalogue = allcatalogues[int(catalogue_number)]
	X, Y, e1_temp, e2_temp = np.loadtxt(catalogue, unpack=True, usecols=(0,1,5,6))
	# You have actually read in ra,dec above which both run from [-3,+3]. Move up/across to origin.
	X += abs(np.min(X))
	Y += abs(np.min(Y))

	
else: # it's the 100 sqdeg run.
	print("Reading in the %s sqdeg shear catalogue for LOS %s" %(sqdeg, los))
	# Pxls already in mask frame
	X, Y, e1_temp, e2_temp = np.loadtxt('%s/Mass_Recon/%s/%s.%sGpAM.LOS%s_Xm_Ym_e1_e2_w.dat'%(data_DIR, DIRname, name, gpam, los), unpack=True, usecols=(0, 1, 2,3))






######################## SHAPE NOISE #############################
# Decide if you want a noise-only, shear+noise, or shear-only calculation
if SN == 'ALL' or SN == 'All' or SN == 'all' or 'Cycle' in DIRname:
	
	if float(gpam) == 3.32: 
		SN_level=0.28
	else:
		SN_level=0.29
	print("SN_level is %s" %SN_level)
	# make it so the SN for a given LOS and cosmol is ALWAYS the same, using np.random.seed

	if 'Cycle' in DIRname:
		# get noise realisation number
		intlos = int(los.split('n')[0])	
		ncycle = int(los.split('n')[-1])
		seed1 = intlos + ncycle*100
		seed2 = intlos + 201 + ncycle*100
		print("seed1 and seed2 are %s and %s"%(seed1,seed2))
	else:
		seed1 = int(los)
		seed2 = int(los) + 201 # made it so all cosmol have same SN.
								# So only diff in signal is due to cosmol.
        
    #np.random.seed(seed1)
	#np.random.normal(0., SN_level, len(e1_temp)) # Add shape noise
	#np.random.seed(seed2)
	#np.random.normal(0., SN_level, len(e1_temp))

	# 28/05/2019 - edit to make sure ellipticities bounded by -1 and 1
	e1_rng = Generate_Unitary_Shape_Noise(seed1, SN_level)
	e2_rng = Generate_Unitary_Shape_Noise(seed1, SN_level)

	if SN == 'ALL' or SN == 'All' or SN == 'all':
		e1_temp = 0. # noise-only
		e2_temp = 0.

elif float(SN) == 0.:
	e1_rng = 0. # shear-only
	e2_rng = 0.

else:
	# make it so the SN for a given LOS is ALWAYS the same, using np.random.seed
	#np.random.seed(int(los))
	#e1_rng = np.random.normal(0., float(SN), len(e1_temp)) # Add shape noise
	#np.random.seed(int(los)+201)
	#e2_rng = np.random.normal(0., float(SN), len(e1_temp))

	# 28/05/2019 - edit to make sure ellipticities bounded by -1 and 1
	e1_rng = Generate_Unitary_Shape_Noise(int(los), float(SN))
	e2_rng = Generate_Unitary_Shape_Noise(int(los)+201, float(SN))		


######################## SHAPE NOISE #############################

# 28/05/2019 - edit to do more accurate contribution of SN to e_obs
e_rng = e1_rng + 1j*e2_rng
e_temp = e1_temp + 1j*e2_temp
e_obs = (e_rng + e_temp) / (1+ e_rng*np.conj(e_temp))
e1 = np.real(e_obs)
e2 = np.imag(e_obs)

#e1 = e1_temp + e1_rng
#e2 = e2_temp + e2_rng


# Convert coords to mask frame
X = (X*float(PS))/PSm
Y = (Y*float(PS))/PSm




# Now Apply Mask if name specifies to.
if mask == 'mask' or mask=='G9mask' or mask=='W3mask':
	 X, Y, e1, e2 = Mask_Shear(X, Y, e1, e2, mask_filename) 

	
No_gals = len(e1)
Weight = np.zeros(No_gals)+1. # Give all gals a weight of unity
np.savetxt('%s/Mass_Recon/%s/%s.%sGpAM.LOS%s_Xm_Ym_e1_e2_w.dat'%(data_DIR, DIRname, name, gpam, los), np.c_[X, Y, e1, e2, Weight], header = '%s 5'%(No_gals), comments='')


# FINALY...
# (03/02/2017) 
# Save another catalogue with phoney ra,dec instead of x,y
# This is so we can load the mock catalogues as healpix maps, and see if 
# healpix VS Ludo mass reconstruction give the same results.
if sqdeg == 100 or sqdeg == 60:
	ra = X*PSm
	dec = -0.5*(new_sizeY*PSm) + Y*PSm # between (-5,+5) or (-3.95,+3.95)
									   # depending on if sqdeg is 100 or 60 deg^2
	
	np.savetxt('%s/Mass_Recon/%s/%s.%sGpAM.LOS%s_ra_dec_e1_e2_w.dat'%(data_DIR, DIRname, name, gpam, los), np.c_[ra, dec, e1, e2, Weight], comments='')






