import numpy as np
from astropy.io import fits
import time
import os
from scipy.stats import binned_statistic_2d
from ClassWarfare import Filter_Input, Format_2Darray
# These classes are used by the Kappa2Shear function


def MeanQ_VS_XY(Q, w, m, X,Y,num_XY_bins):
        # we want the weighted mean of Q and also calibrated.
        # Calculate the sum 2D binned value of Q*w and m*w, and then divide
        sumQw_grid, yedges, xedges, binnum = binned_statistic_2d(Y, X, Q*w, statistic='sum', bins=num_XY_bins)
        sum_mw_grid, yedges, xedges, binnum = binned_statistic_2d(Y, X, m*w, statistic='sum', bins=num_XY_bins)
        AvQ_grid=sumQw_grid/sum_mw_grid

        # Correct bad pixels
        bad_pxls = np.where( np.isfinite(AvQ_grid) == False )
        AvQ_grid[bad_pxls[0], bad_pxls[1]] = 0.
        return AvQ_grid,yedges,xedges,bad_pxls
        


def interpolate2D(X, Y, grid): #(It's linear)
	Xi = np.int64(X) 
	Yi = np.int64(Y) # these round down to integer

	# if any pxl exceeds grid, or is at the egde of grid
        # shift it in 2 pxls, so we can interpolate from this to boundary
	Xi[Xi >= grid.shape[1]-1] = grid.shape[1]-2 # (-2 so Xi+1 below runs okay).
	Yi[Yi >= grid.shape[0]-1] = grid.shape[0]-2

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
 	# This func correctly assembles DIRname irrespective of which cat (hi/lo) is inputted first
 	# and otherwise returns data corresponding to the 1st catalogue (A)
 	variableA = Filter_Input(input_args[:-1])                     # omitting the 2nd paramfile
 	variableB = Filter_Input(input_args[0:2]+[input_args[-1]]+input_args[3:5])

 	variableA.Filter()
 	variableB.Filter()

 	name, gpam, DIRnameA, SS, sigma, SN, mask, z, PS, sqdeg, zlo_A, zhi_A, ThBins, OATH, los, los_end = variableA.Unpack_Sims()
 	_,_,DIRnameB,_,_,SN2,_,_,_,_, zlo_B,zhi_B,_,_,_,_ = variableB.Unpack_Sims()

 	# check which DIRname corresponds to the lower redshift bin.
 	# need this info to correctly assemble the DIRname, so the lower zbin-info is always first.
 	if float(zlo_A) < float(zlo_B):  # Low-z paramfile is first
 		DIRname1,DIRname2 = DIRnameA,DIRnameB
 		zlo,zhi = zlo_A,zhi_A
 	else:                            # High-z paramfile is first
 		DIRname1,DIRname2 = DIRnameB,DIRnameA
 		zlo,zhi = zlo_B,zhi_B

 	DIRname = DIRname1.split('ZBcut')[0] + 'ZBcut%s-%s_X_ZBcut' %(zlo,zhi) + DIRname2.split('ZBcut')[-1]
 	return name, gpam, DIRname, SS, sigma, SN, mask, z, PS, sqdeg, zlo_A, zhi_A, ThBins, OATH, los, los_end


def Shuffle_filenames_LOSnR( Shuffle_Config, filenames, R ):
        # Read in the matrix which describes how the LOS & R are shuffled;
        # this depends on the ID used to shuffle: Shuffle_Config (fiducially set to 0):
        parent_dir = '/home/bengib/Clipping_Pipeline/Mass_Recon/Shuffled_SLICS-K1000-Mosaic_Configs'

        # Determine if there's multiple noise realns:
        if 'Cycle' in filenames[0]:
                # figure out how many noise realns there are:
                Nnoise = len( [i for i in filenames if 'LOS292' in i] ) # how many LOS292 exist for this R
                shffld_lr = np.load('%s/Shuffled_Config%s_Nnoise%s.npy' %(parent_dir,Shuffle_Config,Nnoise))
                NLOS = len( [i for i in filenames if 'R%sn0'%R in i] )  # how many LOS for this R and n=0
        else:
                shffld_lr = np.load('%s/Shuffled_Config%s.npy' %(parent_dir,Shuffle_Config))
                NLOS = len(filenames)
                Nnoise = 1
                
        # Re-shuffle the filenames in accordance with the (Nnoise x) 217x18 shuffle matrix
        files_shffld = []
        count = 0
        for i in range(NLOS):
                for j in range(Nnoise): 
                        los_num = filenames[count].split('LOS')[-1].split('R')[0]
                        # replace the LOS num with another from the matrix
                        # (which is different for each of the 18 regions):
                        if 'Cycle' in filenames[0]:
                                # Could use j (noise idx) to pick which layer of shuffle matrix to get LOS from
                                # but this assumes filenames really are ordered by noise real. num.
                                # safer to use nn from filename itself.
                                nn = int( filenames[count].split('LOS%sR%sn' %(los_num,R) )[-1].split('.')[0] )
                                tmp_file = filenames[count].replace('LOS%sR%sn%s' %(los_num,R,nn),
                                                                    'LOS%sR%sn%s'%(int(shffld_lr[nn,i,R-1]),R,nn) )
                        else:
                                tmp_file = filenames[count].replace('LOS%sR%s' %(los_num,R),
                                                                    'LOS%sR%s'%(int(shffld_lr[i,R-1]),R) )

                        files_shffld.append( tmp_file )
                        count+=1

        return files_shffld

