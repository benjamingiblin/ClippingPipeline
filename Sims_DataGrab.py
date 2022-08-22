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
from FunkShins import interpolate2D, Mask_Shear, Lower_Res_Mask, Combine_zbin_DIRname

if len(sys.argv)-1 > 4 and "param_files" in sys.argv[-1]:
	# 2  paramfiles have given, it's a combination/cross redshift problem (change input/output DIRname)
	name, gpam, DIRname, SS, sigma, SN, mask, z, PS, sqdeg, zlo, zhi, ThBins, OATH, los, los_end = Combine_zbin_DIRname( sys.argv )
else:
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


if sqdeg==36:
	num_rt_cats = 5		# the number of ray tracing catalogues per realisation
	DH10_datadir='/data/bengib/Clipping_SimsLvW//DH10_Mocks/FaLCoNS/'
	# get the cosmol number
	cosmol = int(DIRname.split('_Cosmol')[-1])

	# get the los number + noise cycle number if applicable
	intlos = int(los.split('n')[0])		
	# get the realisation number
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


# Need to make it so the SN for a given LOS and cosmol is ALWAYS the same,
# but different for each redshift bin,                          
# using np.random.seed()
if type(eval(zlo)) == float and type(eval(zhi)) == float:
	zfactor = int( (eval(zlo)+eval(zhi)) *10000 ) # different for each zbin, and not repeated for any LOS.
else:
        zfactor = 0


######################## SHAPE NOISE #############################
# Decide if you want a noise-only, shear+noise, or shear-only calculation
if SN == 'ALL' or SN == 'All' or SN == 'all' or 'Cycle' in DIRname:

        # Set the noise level (not specified in param_file if doing Cycle):
        if 'KiDS1000' in DIRname:
                # Must set SN to values measured in each KiDS1000 bin.
                bin_edges = np.array([0.1, 0.3, 0.5, 0.7, 0.9, 1.2])
                sigma_e_values = [0.270, 0.258, 0.273, 0.254, 0.270]
                
                # Note: this only works for the 5 zbins used for cosmic shear.
                idx_sig = np.where(float(zlo) == bin_edges)[0][0]
                SN_level = sigma_e_values[idx_sig]
        else:
                SN_level=0.28 
                        
        print("SN_level is %s" %SN_level)

        if 'Cycle' in DIRname:
                # get noise realisation number
                intlos = int(los.split('n')[0])	
                ncycle = int(los.split('n')[-1])
                seed1 = intlos + ncycle*100 + zfactor
                seed2 = intlos + 2001 + ncycle*100 + zfactor
                print("seed1 and seed2 are %s and %s"%(seed1,seed2))
        else:
                seed1 = int(los) + zfactor + 19
                seed2 = int(los) + 2001 + zfactor + 21
                                                  # made it so all cosmol have same SN.
		                                  # So only diff in signal is due to cosmol.
                                                  # Also noise maps have different seed to SLICS.
                print("seed1 and seed2 are %s and %s"%(seed1,seed2))
                e1_temp = np.zeros(len(X)) # noise-only 
                e2_temp = np.zeros(len(X))
                
        
        # 28/05/2019 - edit to make sure ellipticities bounded by -1 and 1
        e1_rng = Generate_Unitary_Shape_Noise(seed1, SN_level)
        e2_rng = Generate_Unitary_Shape_Noise(seed2, SN_level)


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
        seed1 = int(los) + zfactor
        seed2 = int(los) + 2001 + zfactor
        e1_rng = Generate_Unitary_Shape_Noise(seed1, float(SN))
        e2_rng = Generate_Unitary_Shape_Noise(seed2, float(SN))		
        print("seed1 and seed2 are %s and %s"%(seed1,seed2))

        
######################## SHAPE NOISE & IAs #############################
# 28/05/2019 - edit to do more accurate contribution of SN to e_obs:
# e = e1 + j*e2
e_rng = e1_rng + 1j*e2_rng     # noise
e_temp = e1_temp + 1j*e2_temp  # shear

# complex addition of shear and noise:
e_obs = (e_rng + e_temp) / (1+ e_rng*np.conj(e_temp))

#e1 = e1_temp + e1_rng
#e2 = e2_temp + e2_rng

if "IA" in DIRname:
	IA_amp = float(DIRname.split('IA')[-1].split('_')[0])
	IA1, IA2 = np.loadtxt('%s/Mass_Recon/%s/%s.%sGpAM.LOS%s_IA1_IA2.dat'%(data_DIR, DIRname, name, gpam, los),
                              usecols=(0,1), unpack=True)
	e_IA = IA1*IA_amp + 1j*IA2*IA_amp # complex IA
       
	# complex addition of e_obs and e_IA
	e_obs = (e_IA + e_obs) / (1+ e_IA*np.conj(e_obs))

e1 = np.real(e_obs)
e2 = np.imag(e_obs)

# Convert X,Y (in arcmin) to coords to mask frame (PSm is in deg/pxl)
X = X/(PSm*60.)
Y = Y/(PSm*60.)


# Now Apply Mask if name specifies to.
if mask == 'mask' or mask=='G9mask' or mask=='W3mask':
	 X, Y, e1, e2 = Mask_Shear(X, Y, e1, e2, mask_filename) 

	
No_gals = len(e1)
Weight = np.zeros(No_gals)+1. # Give all gals a weight of unity
np.savetxt('%s/Mass_Recon/%s/%s.%sGpAM.LOS%s_Xm_Ym_e1_e2_w.dat'%(data_DIR, DIRname, name, gpam, los), np.c_[X, Y, e1, e2, Weight], header = '%s 5'%(No_gals), comments='')








