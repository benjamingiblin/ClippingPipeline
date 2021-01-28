# 24/08/16 B. M. Giblin, PhD student, Edinburgh
# Read in the "residual"-kappa map. Convert it to a 
# "residual"-shear map.



import sys
import numpy as np
from astropy.io import fits
import random
from subprocess import call
import os
import time
from os import getcwd

pipeline_DIR='/home/bengib/Clipping_SimsLvW/'
data_DIR='/data/bengib/Clipping_SimsLvW/'

classdir = pipeline_DIR + "/ShowSumClass"
sys.path.insert(0, classdir) # add directory in which classes & functions 
							 # are defined to the python path
from ClassWarfare import Filter_Input, Format_2Darray
from FunkShins import interpolate2D, Save_FITS, Kappa2Shear



variable = Filter_Input(sys.argv)
variable.Filter()
RUN = sys.argv[1]

if RUN == 'Sims_Run':
	name, gpam, DIRname, SS, sigma, SN, mask, z, PS, sqdeg, zlo, zhi, ThBins, OATH, los, los_end = variable.Unpack_Sims()
	combined_name = '%s.%sGpAM.LOS%s.SS%s'%(name,gpam,los,SS)
	clipfilepath = '%s/Clipping_K/Clip_Thresholds/ClipThreshold_%ssigma.%sSqdeg.SS%s.txt'%(pipeline_DIR,sigma,sqdeg,SS)
	unclipped_file = '%s/Mass_Recon/%s/%s.%sGpAM.LOS%s_Xm_Ym_e1_e2_w.dat'%(data_DIR,DIRname,name,gpam,los)
else:
	DIRname, Blind, SS, sigma, zlo, zhi, ThBins, OATH, Field = variable.Unpack_KiDS()
	combined_name = '%s.Blind%s.SS%s'%(Field,Blind,SS)
	clipfilepath = '%s/Clipping_K/Clip_Thresholds/ClipThreshold_%ssigma.Blind%s.SS%s.txt'%(pipeline_DIR, sigma, Blind, SS)
	unclipped_file = '%s/Mass_Recon/%s/%s.Blind%s_Xm_Ym_e1_e2_w.dat'%(data_DIR,DIRname,Field,Blind)


indir = '%s/Mass_Recon/%s/'%(data_DIR, DIRname)
outdir = '%s/Clipping_K/%s/'%(data_DIR, DIRname)




t0 = time.time()

############################## STEP 1: CLIP THE KAPPA MAP #############################


	
wholename = '%s%s.Ekappa.fits'%(indir,combined_name)

print("opening Fits files")
kappa = fits.open('%s'%wholename)
kappa_c = fits.open('%s'%wholename) #clipped 
kappa_d = fits.open('%s'%wholename) #delta K




# Check if kappa clip threshold already exists in save file
# If not, calculate it and save
if list(str(sigma))[0] =='X':
	f = open('%s/Clipping_K/Clip_Thresholds/%s_Threshold' %(pipeline_DIR,sigma))
	clip = f.readline()
	clip = float(clip)

else:
	# Check if previous clip theshold exists
	# clipfilepath defined (above) differently for KiDS and Sims.
	if os.path.isfile('%s'%clipfilepath):
		f = open('%s'%clipfilepath)
		clip = f.readline()
		clip = float(clip)
	else:
		kav = np.mean(kappa[0].data) 
		kstd = np.std(kappa[0].data)
		clip = kav + float(sigma)*kstd
		f = open('%s'%clipfilepath, 'w')
		f.write('%s \n'%clip)



print("Got the clipping threshold. It is %s. Now doing the clipping." %clip)
# CLIPPING
yc, xc = np.where( kappa[0].data> clip) # The pxls > clip
kappa_c[0].data[yc,xc] = clip
counter = len(yc) # The number of clipped pxls




# Record the fraction of pxls that were clipped
f = open('%s/Correlation_Function/%s/Clipped_PxlFrac.%s.%ssigma.txt'%(pipeline_DIR, DIRname,combined_name, sigma), 'w')

f.write('The fraction of pxls clipped in the kappa map = \n')
frac = float(counter)/len(np.where(kappa[0].data != 0.)[0] ) # Doesn't include masked regions in calculation
#print('Fraction of pixels clipped is %s'%frac)
f.write('%s \n'%frac)
f.write('When clipped at %s sigma with smoothing scale %s pixels\n' %(sigma,SS) )
f.close()
	

# Residual kappa
kappa_d[0].data = kappa[0].data - kappa_c[0].data

# Save clipped and residual kappa maps
#kappa_d.writeto('%s%s.deltak.CLIP_%ssigma.fits'%(outdir,combined_name,sigma), output_verify='ignore', clobber=True)
#kappa_c.writeto('%s%s.ClippedK.CLIP_%ssigma.fits'%(outdir,combined_name,sigma), output_verify='ignore', clobber=True)


# Record the VOLUME of the clipped peaks - NOTE, PXL scale not included in calc. Include this by hand if you want to compare different pxl scales.
f = open('%s/Correlation_Function/%s/Clipped_PxlVol.%s.%ssigma.txt'%(pipeline_DIR, DIRname,combined_name, sigma), 'w')
f.write('The Volme of the peaks clipped in the kappa map = \n')
f.write('%s \n'%np.sum(kappa_d[0].data))
f.write('When clipped at %s sigma with smoothing scale %s pixels\n' %(sigma,SS) )
f.close()
	

kappa.close()
kappa_c.close()




t1 = time.time()
print("time taken for clipping is %.2f seconds" %(t1-t0))










############################## STEP 2: CONVERT DELTA_KAPPA TO DELTA_SHEAR #############################

# In order to avoid ringing effects with FT, pad kappa array with zeros. Run FT
# and cut the padding regions out at the end. 
# Pad up to 8192 --> smallest power of 2 larger than images.
nX = 8192
nY = 8192

delta_e1_grid, delta_e2_grid = Kappa2Shear(kappa_d[0].data, '%s/Clipping_K/Filter1_Kappa2Shear.fits'%pipeline_DIR, '%s/Clipping_K/Filter2_Kappa2Shear.fits'%pipeline_DIR, nX, nY)
#Save_FITS(delta_e1_grid, '%s%s.deltashear1_Re.CLIP_%ssigma.fits'%(outdir,combined_name,sigma))
#Save_FITS(delta_e2_grid, '%s%s.deltashear2_Re.CLIP_%ssigma.fits'%(outdir,combined_name,sigma))

t2 = time.time()
print("time taken for converting delta-kappa to delta-shear is %.2f seconds" %(t2-t1))






############################## STEP 3: OBTAIN CLIPPED SHEAR #############################



# shear_clipped = shear_orig - shear_resid 
t2a = time.time()
Xm, Ym, e1, e2, w = np.loadtxt(unclipped_file, skiprows=1, unpack=True) # skip line with <np. gals> <no. cols>
t2b = time.time()
print("Reading in the orig shear values took %.2f seconds" %(t2b - t2a))


# Just so interpolation doesn't break... append grids with final rows and columns
# and move maximum pxl value down 1/100th of a pxl:
Xm[Xm==Xm.max()] = Xm.max()-0.01
Ym[Ym==Ym.max()] = Ym.max()-0.01

delta_e1_grid = np.c_[ delta_e1_grid, delta_e1_grid[:,-1] ] 
delta_e1_grid = np.r_[ delta_e1_grid, [delta_e1_grid[-1,:]] ] 

delta_e2_grid = np.c_[ delta_e2_grid, delta_e2_grid[:,-1] ] 
delta_e2_grid = np.r_[ delta_e2_grid, [delta_e2_grid[-1,:]] ] 

# the delta-shear at each of the coordinates at which unclipped is defined
delta_e1_Orig = interpolate2D(Xm, Ym, delta_e1_grid)
delta_e2_Orig = interpolate2D(Xm, Ym, delta_e2_grid)

e1_clipped = e1 - delta_e1_Orig
e2_clipped = e2 - delta_e2_Orig

t2c = time.time()
print("Saving the clipped shear values")
np.savetxt('%s/%s.rCLIP_%ssigma.Xm_Ym_e1c_e2c.asc'%(outdir,combined_name,sigma), np.c_[Xm, Ym, e1_clipped, e2_clipped, w])
t2d = time.time()
print("Saving the clipped shear values took %.2f seconds" %(t2d - t2c))


SaveDelta = "N"
if SaveDelta == "Y":
	np.savetxt('%s/%s.deltaCLIP_%ssigma.Xm_Ym_e1c_e2c.asc'%(outdir,combined_name,sigma), np.c_[Xm, Ym, delta_e1_Orig, delta_e2_Orig, w])	


t3 = time.time()
print("time taken for calculating clipped shear is %.2f seconds" %(t3-t2))



kappa_d.close()

print("TOTAL TIME TAKEN FOR GET_CLIPPED_SHEAR.py IS %.2f seconds" %(t3 - t0))






















