import numpy as np
from astropy.io import fits
import time
import os
from ClassWarfare import Filter_Input, Format_2Darray
# These classes are used by the Kappa2Shear function





def interpolate2D(X, Y, grid): #(It's linear)
	Xi = X.astype(np.int)
	Yi = Y.astype(np.int) # these round down to integer

	VAL_XYlo = grid[Yi, Xi] + (X - Xi)*( grid[Yi, Xi+1] - grid[Yi, Xi] )
	VAL_XYhi = grid[Yi+1,Xi] + (X - Xi)*( grid[Yi+1,Xi+1] - grid[Yi+1, Xi] )
	VAL_XY = VAL_XYlo + (Y - Yi)*( VAL_XYhi - VAL_XYlo )  
	
	return VAL_XY




def Mask_Shear(X, Y, e1, e2, mask_filename):
	Xi = X.astype(np.int)
	Yi = Y.astype(np.int) # these round down to integer

	Mask = fits.open(mask_filename)

	# avoid out of bounds pxls
	if Yi.max() == Mask[0].data.shape[0]:
		Yi[ np.where( Yi == Yi.max() )[0] ] = Yi.max()-1

	if Xi.max() == Mask[0].data.shape[1]:
		Xi[ np.where( Xi == Xi.max() )[0] ] = Xi.max()-1

	MaskVals_4Extractions = Mask[0].data[Yi, Xi]
	Usable_indi = np.where(MaskVals_4Extractions== 0) 
	return X[Usable_indi], Y[Usable_indi], e1[Usable_indi], e2[Usable_indi]






######################### THE FOLLOWING IS ALL INVOLVED IN CONVERTING KAPPA 2 SHEAR ###############################


def Save_FITS(mapp, filename):
	mapp_map = fits.PrimaryHDU()
	mapp_map.data = mapp
	mapp_map.writeto(filename, output_verify='ignore', clobber=True)




# Takes as arguments, addresses of 2 filters + X and Y arrays defining 
# the dimensions and 0-point of the filter. See Kappa2Shear function
# for default.
def Make_Filters(filter_name1, filter_name2, array1, array2):

	t1 = time.time()
		
	filt1 = np.empty([len(array1), len(array2)])
	filt2 = np.empty([len(array1), len(array2)])
	print("Making kappa-to-shear filters")
	y=0 # --> j
	for j in array1:
		x=0 # --> i
		for i in array2:
			filt1[y, x] = ((i**2. -j**2.)/(i**2. + j**2.)) 
			filt2[y, x] =  (2.*i*j/(i**2. + j**2.)) 
			x+=1
		y+=1		
	t2 = time.time()
	print("Making the filter took %.2f seconds" %(t2 - t1))

	# Save filter
	Save_FITS(filt1, filter_name1)
	Save_FITS(filt2, filter_name2)
	return





# filt1 and 2 are the filename addresses of the filters
# nX and nY are the dimensions of the PADDED UP kappa map
# i.e. should be a power of 2.
def Kappa2Shear(kap, filt1, filt2, nX, nY):

	t1 = time.time()
	print("Padding the kappa array and computing the FFT")
	kappa_padded = Format_2Darray(kap).Pad_2Darray(nX, nY, 0.) # Calls class defined in same directory
	print("Now doing the FFt")
	kappaFFT = np.fft.fft2(kappa_padded)
	kappaFFTshift = np.fft.fftshift(kappaFFT)
	t2 = time.time()
	print("Computing the FFT of kappa took %.2f" %(t2 - t1))



	# Make filters if they don't exist
	if os.path.exists(filt1) != True or os.path.exists(filt2) != True:
		# Always pass y-axis array first.
		Make_Filters(filt1, filt2, np.linspace(-nY/2, nY/2, nY), np.linspace(-nX/2, nX/2, nX))


	# Argument is the address of the filter
	def K2S_Components_Conversion(filter_name):

		# Load the filter
		filt_map = fits.open(filter_name)
		filt = (filt_map[0].data)

		# Dont need -1 here, because of way filters are defined.
		# Note 2 self: tried dot product, deffo doesn't work.
		shearfft = filt*kappaFFTshift
		
		shear = np.fft.ifft2( np.fft.ifftshift(shearfft) ) # shift 0-freqs back to edge and IFFT
	
		# Remove the padded regions --> return to original size.
		shear_unpad = Format_2Darray(shear).Cut_Edges(len(kap[0,:]), len(kap[:,0]))


		filt_map.close()
	
		return shear_unpad.real

	
	shear1 = K2S_Components_Conversion(filt1)
	shear2 = K2S_Components_Conversion(filt2)

	return shear1, shear2


###############################################################################################################


# Sort a NLOS * no_bins array into a groups * no_bins array
# Where each column element is the MEAN of 'groups' column elements 
def Sort_Array_IntoGroups(CF_array, groups):
	NLOS = len(CF_array[:,0])
	no_bins = len(CF_array[0,:]) 

	# sort LOS into groups of 4
	if NLOS % groups != 0:
		new_length = int(NLOS/groups) + 1
	else:
		new_length = int(NLOS/groups)

	Sorted_CF_array = np.empty([ new_length, no_bins ])
	
	for i in range(0, no_bins):
		for j in range(0, int(NLOS/groups)):
			Sorted_CF_array[j,i] = np.mean( CF_array[groups*j:groups*(j+1),i] )
		if NLOS % groups != 0:
			Sorted_CF_array[-1,i] = np.mean( CF_array[-1*(NLOS%groups):,i] )

	return Sorted_CF_array



# Read in the name of a mask + its resolution in arcmin, calculate the area and the unmasked area
def Calc_Unmasked_Area(mask_name, res_in_arcmin):    
    mask = fits.open(mask_name)
    No_Pxls = (len(mask[0].data[0,:])*len(mask[0].data[:,0]))
    Tot_Area = No_Pxls * (float(res_in_arcmin)/60.)**2.
    Unmasked_Frac = 1. - float(len(np.where(mask[0].data > 0)[0])) / No_Pxls
    return Unmasked_Frac, Tot_Area






# Read in a 2D Mask array and lower its resolution
def Lower_Res_Mask(Mask_HiRes, new_sizeX, new_sizeY):

	pxls = np.where(Mask_HiRes != 0)
	ym = pxls[0]
	xm = pxls[1]

	Mask_HiRes[ym, xm]=1 # change all the masked pxls to have value of 1


	# Now lower the resolution of this thing to new_sizeX*new_sizeY

	Mask_LoRes = fits.PrimaryHDU() # Make a new fits image object
	Mask_LoRes.data = np.empty([new_sizeY, new_sizeX])
	
	# Averaging the pixels in blocks of jump^2 pxls^2 
	jump = int(len(Mask_HiRes[0,:])/len(Mask_LoRes.data[0,:]) ) # =28
	

	# For loop goes through pxls of lo res mask: new_sizeX*new_sizeY
	for x in range(0, len(Mask_LoRes.data[0,:]) ):	
		for y in range(0, len(Mask_LoRes.data[:,0]) ):
			#print('\n \n Summing x elems, %d to %d \n'%(x*jump, (x+1)*jump))
			#print('Summing y elems, %d to %d \n'%(y*jump, (y+1)*jump))
			avg_pxls = float(round( np.sum( Mask_HiRes[ y*jump:(y+1)*jump, x*jump:(x+1)*jump ] )/(jump*jump) )) # Make this an integer lad....
			#print(avg_pxls)
			Mask_LoRes.data[y,x] = avg_pxls


	masked = np.where(Mask_LoRes.data[:,:] != 0)
	frac = float(len(masked[0]))/(new_sizeX*new_sizeY)
	print('Masked fraction of image is %s'%frac)
	# This prints out ~27%

	return Mask_LoRes.data, masked 









def ReBin_CF(theta, CF, npairs, new_bins):

	# Re-bin a long CF into a smaller no. of bins.
	# IF len(CF) does not divide perfectly into new_bins, this function works best when len(CF) >> new_bins
	# Otherwise you get a lot more contributions to new_CF[0] than the other new_bins

	theta_new = np.zeros(new_bins)
	CF_new = np.zeros(new_bins) 
	npairs_new = np.zeros(new_bins) # the number of gal pairs that go into each new bin.

	collect = len(CF) / new_bins # the number of original bins that go into each of the new bins
	r = len(CF) % new_bins 
	start_collect = collect + r # the number of elems at the start that will be sorted into 
												  # the first bin ... this is to account for if len(CF) is not divisible by new_bins
	
	
	# sort the first 'start_collect' no. of elems into bin 0 of the new CF
	# It looks like Athena is returning the 10 ** mean-of-log-theta
	# So for consistency, we should return: 10 ** weight-mean-of-log-theta
	CF_new[0] = np.sum(npairs[0:start_collect] * CF[0:start_collect]) / np.sum(npairs[0:start_collect])
	theta_new[0] = 10 ** (np.sum(npairs[0:start_collect] * np.log10(theta[0:start_collect]) ) / np.sum(npairs[0:start_collect]))
	npairs_new[0] = np.sum(npairs[0:start_collect])

	for i in range( 1, new_bins ):
		CF_new[i] = np.sum(npairs[i*collect+r : (i+1)*collect+r] * CF[i*collect+r : (i+1)*collect+r]) / np.sum(npairs[i*collect+r : (i+1)*collect+r])
		theta_new[i] = 10 ** (np.sum(npairs[i*collect+r : (i+1)*collect+r] * np.log10(theta[i*collect+r : (i+1)*collect+r]) ) / np.sum(npairs[i*collect+r : (i+1)*collect+r]))
		npairs_new[i] = np.sum(npairs[i*collect+r : (i+1)*collect+r])

	return theta_new, CF_new, npairs_new



def Combine_zbin_DIRname(input_args):

	variable1 = Filter_Input(input_args[:-1])                     # omitting the 2nd paramfile
	variable2 = Filter_Input(input_args[0:2]+[input_args[-1]]+input_args[3:5])

	variable1.Filter()
	variable2.Filter()

	name, gpam, DIRname1, SS, sigma, SN, mask, z, PS, sqdeg, zlo, zhi, ThBins, OATH, los, los_end = variable1.Unpack_Sims()
	_,_,DIRname2,_,_,SN2,_,_,_,_, zlo2,zhi2,_,_,_,_ = variable2.Unpack_Sims()

	DIRname = DIRname1.split('ZBcut')[0] + 'ZBcut%s-%s_X_ZBcut' %(zlo,zhi) + DIRname2.split('ZBcut')[-1]
	return name, gpam, DIRname, SS, sigma, SN, mask, z, PS, sqdeg, zlo, zhi, ThBins, OATH, los, los_end




