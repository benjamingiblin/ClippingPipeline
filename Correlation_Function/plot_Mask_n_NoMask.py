# 05/03/2017, B. M. Giblin, PhD Student Edinburgh
# Compare masked and unmasked CFs. Type in the param file for the masked measurement
# The unmasked will be determined and read in automatically.




import pylab as plt
import numpy as np
import sys
import os
import matplotlib
from subprocess import call
import glob
from os import getcwd
from matplotlib import rc
from scipy.interpolate import interp1d

overall_DIR = getcwd()
classdir = getcwd() + "/ShowSumClass"
sys.path.insert(0, classdir) # add directory in which classes & functions 
							 # are defined to the python path
from ClassWarfare import Filter_Input, Handle_CF_Files
from FunkShins import Sort_Array_IntoGroups, Calc_Unmasked_Area

variable = Filter_Input(sys.argv)
variable.Filter()
RUN = sys.argv[1]




rc('text',usetex=True)
rc('font',size=18)
rc('legend',**{'fontsize':18})
rc('font',**{'family':'serif','serif':['Computer Modern']})







if RUN == 'Sims_Run':

	# Read in directory which has 'Mask' specified
	name, gpam, DIRname, SS, sigma, SN, mask, z, PS, sqdeg, zlo, zhi, ThBins, OATH, los_start, los_end = variable.Unpack_Sims()
	
	# additionally read in NoMask data...
	DIRname2 = DIRname.split('W3Mask')[0] + 'NoMask' + DIRname.split('W3Mask')[-1]
	name2 = name.split('W3Mask')[0] + 'test'

	# ... AND the W3Mask data
	DIRname3 = DIRname.split('Mask')[0] + 'W3Mask' + DIRname.split('Mask')[-1]
	name3 = name.split('Mask')[0] + 'W3Mask'


	# ... AND the G9 & unmasked NOISE-ONLY calculation
	DIRname4 = '100Sqdeg_NOISE_W3Mask_8.53GpAM_zKiDS_ZBcutNone'
	name4 = 'NOISE_W3Mask'

	DIRname5 = '100Sqdeg_NOISE_NoMask_8.53GpAM_zKiDS_ZBcutNone'
	name5 = 'NOISE_test'


	indir = '%s/Correlation_Function/%s/ThBins%s' %(overall_DIR, DIRname, ThBins)
	indir2 = '%s/Correlation_Function/%s/ThBins%s' %(overall_DIR, DIRname2, ThBins)
	indir3 = '%s/Correlation_Function/%s/ThBins%s' %(overall_DIR, DIRname3, ThBins)	
	indir4 = '%s/Correlation_Function/%s/ThBins%s' %(overall_DIR, DIRname4, ThBins)	
	indir5 = '%s/Correlation_Function/%s/ThBins%s' %(overall_DIR, DIRname5, ThBins)
	
	# G9 masked
	files_uc = glob.glob('%s/%s.%sGpAM.LOS*.ORIG.CorrFun.asc'%(indir,name,gpam))
	files_c = glob.glob('%s/%s.%sGpAM.LOS*.SS%s.rCLIP_%ssigma.CorrFun.asc'%(indir,name,gpam,SS,sigma))
	NLOS=len(files_uc)

	# unmasked
	files_uc2 = glob.glob('%s/%s.%sGpAM.LOS*.ORIG.CorrFun.asc'%(indir2,name2,gpam))
	files_c2 = glob.glob('%s/%s.%sGpAM.LOS*.SS%s.rCLIP_%ssigma.CorrFun.asc'%(indir2,name2,gpam,SS,sigma))

	# W3 masked
	#files_uc3 = glob.glob('%s/%s.%sGpAM.LOS*.ORIG.CorrFun.asc'%(indir3,name3,gpam))
	#files_c3 = glob.glob('%s/%s.%sGpAM.LOS*.SS%s.rCLIP_%ssigma.CorrFun.asc'%(indir3,name3,gpam,SS,sigma))	


	# NOISE stuff
	files_uc4 = glob.glob('%s/%s.%sGpAM.LOS*.ORIG.CorrFun.asc'%(indir4,name4,gpam))
	files_c4 = glob.glob('%s/%s.%sGpAM.LOS*.SS%s.rCLIP_%ssigma.CorrFun.asc'%(indir4,name4,gpam,SS,sigma))

	files_uc5 = glob.glob('%s/%s.%sGpAM.LOS*.ORIG.CorrFun.asc'%(indir5,name5,gpam))
	files_c5 = glob.glob('%s/%s.%sGpAM.LOS*.SS%s.rCLIP_%ssigma.CorrFun.asc'%(indir5,name5,gpam,SS,sigma))


	outdir='%s/Correlation_Function/%s/ThBins%s/NLOS%s/'%(overall_DIR, DIRname, ThBins, NLOS)
	outdir2='%s/Correlation_Function/%s/ThBins%s/NLOS%s/'%(overall_DIR, DIRname2, ThBins, NLOS)
	outdir3='%s/Correlation_Function/%s/ThBins%s/NLOS%s/'%(overall_DIR, DIRname3, ThBins, NLOS)	
	outdir4='%s/Correlation_Function/%s/ThBins%s/NLOS%s/'%(overall_DIR, DIRname4, ThBins, NLOS)	
	outdir5='%s/Correlation_Function/%s/ThBins%s/NLOS%s/'%(overall_DIR, DIRname5, ThBins, NLOS)	

	combined_name = '%s%s.%sGpAM.NLOS%.0f'%(outdir,name,gpam,NLOS)
	combined_name2 = '%s%s.%sGpAM.NLOS%.0f'%(outdir2,name2,gpam,NLOS)
	combined_name3 = '%s%s.%sGpAM.NLOS%.0f'%(outdir3,name3,gpam,NLOS)
	combined_name4 = '%s%s.%sGpAM.NLOS%.0f'%(outdir4,name4,gpam,NLOS)
	combined_name5 = '%s%s.%sGpAM.NLOS%.0f'%(outdir5,name5,gpam,NLOS)



else:
	
	DIRname, Blind, SS, sigma, zlo, zhi, ThBins, OATH, Dummy = variable.Unpack_KiDS()
	# check if multiple fields are given as input
	arguments = sys.argv[3:]
	Field = []
	for f in arguments:
		if list(f)[0] == 'G':
			Field.append(f)

	indir='%s/Tree_Correlation_Function/%s/ThBins%s/' %(overall_DIR, DIRname, ThBins)
	indir2 = indir # For KiDS, cannot run unmasked run of data - so just load in the masked data again
	indir4='%s/Tree_Correlation_Function/%s_NOISE/ThBins%s/' %(overall_DIR, DIRname, ThBins)
	indir5='%s/Tree_Correlation_Function/%s_NOISE_NoMask/ThBins%s/' %(overall_DIR, DIRname, ThBins) # For KiDS, unmasked noise computed with TreeCorr


	# UNCLIPPED FILES 
	files_uc=[] 	# Mask 
	files_uc2=[]	# For KiDS, cannot run unmasked run of data - so just load in the masked data again
	files_uc4=[]	# Mask and NOISE-ONLY
	files_uc5=[]	# No Mask and NOISE_ONLY
	# CLIPPED FILES	
	files_c=[]	# Mask
	files_c2=[] # For KiDS, cannot run unmasked run of data - so just load in the masked data again

	files_c4=[]	# Mask and NOISE-ONLY
	files_c5=[]	# No Mask and NOISE-ONLY	 


	filename='' # the beginning part of output filenames
				# ='G*G*G*...' depending on number of fields
	print 'You have selected Field(s):'
	for f in Field:
		print '%s '%f
		# UNCLIPPED
		files_uc.append('%s/%s.Blind%s.ORIG.CorrFun.asc' %(indir,f,Blind))
		files_uc2.append('%s/%s.Blind%s.ORIG.CorrFun.asc' %(indir2,f,Blind))
		# CLIPPED
		files_c.append('%s/%s.Blind%s.SS%s.rCLIP_%ssigma.CorrFun.asc' %(indir,f,Blind,SS,sigma))
		files_c2.append('%s/%s.Blind%s.SS%s.rCLIP_%ssigma.CorrFun.asc' %(indir2,f,Blind,SS,sigma))

		# NOISE
		NLOS_NOISE = 25
		for nn in range(0, NLOS_NOISE):
			files_uc4.append('%s/%s_NOISE%s.Blind%s.ORIG.CorrFun.asc' %(indir4,f,nn,Blind))
			files_uc5.append('%s/%s_NOISE%s.Blind%s.ORIG.CorrFun.asc' %(indir5,f,nn,Blind))

			files_c4.append('%s/%s_NOISE%s.Blind%s.SS%s.rCLIP_%ssigma.CorrFun.asc' %(indir4,f,nn,Blind,SS,sigma))	
			files_c5.append('%s/%s_NOISE%s.Blind%s.SS%s.rCLIP_%ssigma.CorrFun.asc' %(indir5,f,nn,Blind,SS,sigma))	
			
		filename = filename + f
	filename_NOISE = filename + '_NOISE%s' %NLOS_NOISE # already has the 'f' added in

	outdir = indir
	outdir2 = indir2
	outdir4 = indir4
	outdir5 = indir5

	combined_name = '%s/%s.Blind%s'%(outdir, filename, Blind)
	combined_name2 = '%s/%s.Blind%s'%(outdir2, filename, Blind)	
	combined_name4 = '%s/%s.Blind%s'%(outdir4, filename_NOISE, Blind)	
	combined_name5 = '%s/%s.Blind%s'%(outdir5, filename_NOISE, Blind)
	
	NLOS  = len(Field)


no_bins = int(ThBins)
def Return_CFpm_And_CFpmArrays(files):

	Class = Handle_CF_Files(files, no_bins)
	if 'Tree' in files[0]:
		col_CF=3
		col_w=8
	else:
		col_CF=1
		col_w=7
	theta, CFp, theta_array, CFp_array, npair_array = Class.Average_CFs(col_CF,col_w)
	theta, CFm, theta_array, CFm_array, npair_array = Class.Average_CFs(col_CF+1,col_w)

	return theta, CFp, CFp_array, CFm, CFm_array, Class 



def Calc_AutoCovMat(CFp, CFp_array, CFm, CFm_array, Class, savename):

	errCFp, CCC_CFp = Class.Calc_Covariance(CFp_array, CFp, CFp_array, CFp, '%s+.'%savename)
	errCFm, CCC_CFm = Class.Calc_Covariance(CFm_array, CFm, CFm_array, CFm, '%s-.'%savename)
	CF= np.vstack(( CFp, CFm ))
	errCF = np.vstack(( np.sqrt(errCFp), np.sqrt(errCFm) ))

	return CF, errCF




####################################### MASKED #########################################
# Unclipped
theta, CFp_uc_tmp1, CFp_uc_array, CFm_uc_tmp1, CFm_uc_array, Class_uc = Return_CFpm_And_CFpmArrays(files_uc) 
#CF_uc, errCF_uc = Calc_AutoCovMat(CFp, CFp_uc_array, CFm, CFm_uc_array, Class, '%s.UCxUC'%(combined_name))
# Clipped
theta, CFp_c_tmp1, CFp_c_array, CFm_c_tmp1, CFm_c_array, Class_c = Return_CFpm_And_CFpmArrays(files_c) 
#CF_c, errCF_c = Calc_AutoCovMat(CFp, CFp_c_array, CFm, CFm_c_array, Class, '%s.SS%s.rCLIP_%ssigma.CxC'%(combined_name,SS,sigma))


###################################### UNMASKED ########################################
# Unclipped
theta, CFp_uc_tmp2, CFp_uc_array2, CFm_uc_tmp2, CFm_uc_array2, Class_uc2 = Return_CFpm_And_CFpmArrays(files_uc2) 
#CF_uc2, errCF_uc2 = Calc_AutoCovMat(CFp, CFp_uc_array2, CFm, CFm_uc_array2, Class, '%s.UCxUC'%(combined_name2))
# Clipped
theta, CFp_c_tmp2, CFp_c_array2, CFm_c_tmp2, CFm_c_array2, Class_c2 = Return_CFpm_And_CFpmArrays(files_c2) 
#CF_c2, errCF_c2 = Calc_AutoCovMat(CFp, CFp_c_array2, CFm, CFm_c_array2, Class, '%s.SS%s.rCLIP_%ssigma.CxC'%(combined_name2,SS,sigma))


###################################### MASKED/UNMASKED COVARIANCE ########################################
###################################### HAS TO BE DONE SEPARATELY TO NOISE REALISATIONS, IN CASE THERES JUST ONE FIELD ########################################
if NLOS == 1: # Covariance estimation will break, so peg some ones onto CF_arrays
	CFp_uc_array = np.vstack(( CFp_uc_array, CFp_uc_array ))
	CFp_c_array = np.vstack(( CFp_c_array, CFp_c_array ))
	CFm_uc_array = np.vstack(( CFm_uc_array, CFm_uc_array ))
	CFm_c_array = np.vstack(( CFm_c_array, CFm_c_array ))

	CFp_uc_array2 = np.vstack(( CFp_uc_array2, CFp_uc_array2 ))
	CFp_c_array2 = np.vstack(( CFp_c_array2, CFp_c_array2 ))
	CFm_uc_array2 = np.vstack(( CFm_uc_array2, CFm_uc_array2 ))
	CFm_c_array2 = np.vstack(( CFm_c_array2, CFm_c_array2 ))

CF_uc, errCF_uc = Calc_AutoCovMat(CFp_uc_tmp1, CFp_uc_array, CFm_uc_tmp1, CFm_uc_array, Class_uc, '%s.UCxUC'%(combined_name))
CF_c, errCF_c = Calc_AutoCovMat(CFp_c_tmp1, CFp_c_array, CFm_c_tmp1, CFm_c_array, Class_c, '%s.SS%s.rCLIP_%ssigma.CxC'%(combined_name,SS,sigma))
CF_uc2, errCF_uc2 = Calc_AutoCovMat(CFp_uc_tmp2, CFp_uc_array2, CFm_uc_tmp2, CFm_uc_array2, Class_uc2, '%s.UCxUC'%(combined_name2))
CF_c2, errCF_c2 = Calc_AutoCovMat(CFp_c_tmp2, CFp_c_array2, CFm_c_tmp2, CFm_c_array2, Class_c2, '%s.SS%s.rCLIP_%ssigma.CxC'%(combined_name2,SS,sigma))





################################# X-COVARIANCE BETWEEN MASKED AND UNMASKED #########################################
# Class name now becomes unimportant
Class = Class_uc
# Unclipped
errCFp_uc_MUM, CCC_CFp_uc_MUM = Class.Calc_Covariance(CFp_uc_array, CF_uc[0,:], CFp_uc_array2, CF_uc2[0,:], '%s.MaskedUCxUnmaskedUC+.'%(combined_name))
errCFm_uc_MUM, CCC_CFm_uc_MUM = Class.Calc_Covariance(CFm_uc_array, CF_uc[1,:], CFm_uc_array2, CF_uc2[1,:], '%s.MaskedUCxUnmaskedUC-.'%(combined_name))
errCF_uc_MUM = np.vstack(( errCFp_uc_MUM, errCFm_uc_MUM ))

# Now for clipped
errCFp_c_MUM, CCC_CFp_c_MUM = Class.Calc_Covariance(CFp_c_array, CF_c[0,:], CFp_c_array2, CF_c2[0,:], '%s.SS%s.rCLIP_%ssigma.MaskedCxUnmaskedC+.'%(combined_name,SS,sigma))
errCFm_c_MUM, CCC_CFm_c_MUM = Class.Calc_Covariance(CFm_c_array, CF_c[1,:], CFm_c_array2, CF_c2[1,:], '%s.SS%s.rCLIP_%ssigma.MaskedCxUnmaskedC-.'%(combined_name,SS,sigma))
errCF_c_MUM = np.vstack(( errCFp_c_MUM, errCFm_c_MUM ))







####################################### NOISE MASKED #########################################
# Unclipped
theta, CFp, CFp_array4, CFm, CFm_array4, Class = Return_CFpm_And_CFpmArrays(files_uc4) 
CF_uc4, errCF_uc4 = Calc_AutoCovMat(CFp, CFp_array4, CFm, CFm_array4, Class, '%s.UCxUC'%(combined_name4))
# Clipped
theta, CFp, CFp_array4, CFm, CFm_array4, Class = Return_CFpm_And_CFpmArrays(files_c4) 
CF_c4, errCF_c4 = Calc_AutoCovMat(CFp, CFp_array4, CFm, CFm_array4, Class, '%s.SS%s.rCLIP_%ssigma.CxC'%(combined_name4,SS,sigma))



###################################### NOISE UNMASKED #######################################
# Unclipped
theta, CFp, CFp_array5, CFm, CFm_array5, Class = Return_CFpm_And_CFpmArrays(files_uc5) 
CF_uc5, errCF_uc5 = Calc_AutoCovMat(CFp, CFp_array5, CFm, CFm_array5, Class, '%s.UCxUC'%(combined_name5))
# Clipped
theta, CFp, CFp_array5, CFm, CFm_array5, Class = Return_CFpm_And_CFpmArrays(files_c5) 
CF_c5, errCF_c5 = Calc_AutoCovMat(CFp, CFp_array5, CFm, CFm_array5, Class, '%s.SS%s.rCLIP_%ssigma.CxC'%(combined_name5,SS,sigma))


###################################### IF KiDS, need to average Noise realisations for each field #######################################
def Av_Noise_Per_KiDS_Field(array):
	temp = np.zeros(no_bins)
	temp_array = np.zeros(no_bins)
	k=0
	for i in range(0, len(Field)):
		for j in range(no_bins):
			temp[j] = np.mean( array[k*NLOS_NOISE:(k+1)*NLOS_NOISE,j] )
		k+=1
		temp_array = np.vstack(( temp_array, temp ))
	temp_array = np.delete(temp_array, (0), axis=0)
	return temp_array

if RUN == 'KiDS_Run':
	# Noise Masked
	CFp_array4 = Av_Noise_Per_KiDS_Field(CFp_array4)
	CFm_array4 = Av_Noise_Per_KiDS_Field(CFm_array4)
	# Noise Unmasked
	CFp_array5 = Av_Noise_Per_KiDS_Field(CFp_array5)
	CFm_array5 = Av_Noise_Per_KiDS_Field(CFm_array5)

		
########################################################################################################################################	






################################# THE CORRECTION, NOISE_MASKED - NOISE_UNMASKED #######################
Correction = CF_c4 - CF_c5
Correctionp_Array = CFp_array4 - CFp_array5
Correctionm_Array = CFm_array4 - CFm_array5

CF_c_Corrected = CF_c - Correction
print CF_c_Corrected
CFp_array_Corrected = CFp_c_array - Correctionp_Array
CFm_array_Corrected = CFm_c_array - Correctionm_Array

# Error on the corrected 
errCFp_c_Corrected, CCC_CFp = Class.Calc_Covariance(CFp_array_Corrected, CF_c_Corrected[0,:], CFp_array_Corrected, CF_c_Corrected[0,:], '%s.SS%s.rCLIP_%ssigma.CorrectedMasked_CxC+.'%(combined_name,SS,sigma))
errCFm_c_Corrected, CCC_CFm = Class.Calc_Covariance(CFm_array_Corrected, CF_c_Corrected[1,:], CFm_array_Corrected, CF_c_Corrected[1,:], '%s.SS%s.rCLIP_%ssigma.CorrectedMasked_CxC-.'%(combined_name,SS,sigma))
errCF_c_Corrected = np.vstack(( np.sqrt(errCFp_c_Corrected), np.sqrt(errCFm_c_Corrected) ))




############################## X-CORRELATION BETWEEN CORRECTED_MASKED_CLIPPED AND UNMASKED ###################################

errCFp_c_MUM_Corrected, CCC_CFp_c_MUM_Corrected = Class.Calc_Covariance(CFp_array_Corrected, CF_c_Corrected[0,:], CFp_c_array2, CF_c2[0,:], '%s.SS%s.rCLIP_%ssigma.CorrectedMaskedCxUnmaskedC+.'%(combined_name,SS,sigma))
errCFm_c_MUM_Corrected, CCC_CFm_c_MUM_Corrected = Class.Calc_Covariance(CFm_array_Corrected, CF_c_Corrected[1,:], CFm_c_array2, CF_c2[1,:], '%s.SS%s.rCLIP_%ssigma.CorrectedMaskedCxUnmaskedC-.'%(combined_name,SS,sigma))

errCF_c_MUM_Corrected = np.vstack(( errCFp_c_MUM_Corrected, errCFm_c_MUM_Corrected ))





# THIS STILL NEEDS WORK
############################## IF Run='KiDS_Run' THEN IT WILL HAVE CALCULATED USELESS ERRORS FROM DATA ###################################
######################################### SCRAP THESE AND LOAD THE ONES FROM THE SIMULATIONS ##############################################
if RUN =='KKiDS_Run':
	ThBins_Sims = 16
	Sims_Filename = '%s/Correlation_Function/100Sqdeg_SN0.29_W3Mask_8.53GpAM_zKiDS_ZBcutNone/ThBins%s/NLOS114/SN0.29_W3Mask.8.53GpAM.NLOS114' %(overall_DIR,ThBins_Sims)
	# Work out unmasked area in Sims to get rescaling right
	Unmasked_Frac, Tot_Area = Calc_Unmasked_Area('%s/W3.16bit.5arcs.reg.Now_100sqdeg.fits'%overall_DIR, 5./60.)
	Rescale = Unmasked_Frac*Tot_Area / 360. # Errors read in are stdev not mean... so no factor of 114 necessary...
	
	# Get the errors then...
	# Signal
	errCFp_uc = np.loadtxt('%s.AverageCF+.asc' %Sims_Filename, usecols=(2,), unpack=True)
	errCFm_uc = np.loadtxt('%s.AverageCF+.asc' %Sims_Filename, usecols=(2,), unpack=True)
	errCF_uc = np.vstack(( errCFp_uc*np.sqrt(Rescale), errCFm_uc*np.sqrt(Rescale) ))
	errCF_uc2 = errCF_uc

	errCFp_c = np.loadtxt('%s.SS112.rCLIP_X1sigma.AverageCF+.asc' %Sims_Filename, usecols=(2,), unpack=True)
	errCFm_c = np.loadtxt('%s.SS112.rCLIP_X1sigma.AverageCF-.asc' %Sims_Filename, usecols=(2,), unpack=True)	
	errCF_c = np.vstack(( errCFp_c*np.sqrt(Rescale), errCFm_c*np.sqrt(Rescale) ))
	errCF_c2 = errCF_c

	
	# The cross-talk between masked and unmasked
	temp = np.load('%s.SS112.rCLIP_X1sigma.CorrectedMaskedCxUnmaskedC+.CovMat.npy'%(Sims_Filename))
	errCFp_c_MUM_Corrected = temp.diagonal()*Rescale
	temp = np.load('%s.SS112.rCLIP_X1sigma.CorrectedMaskedCxUnmaskedC-.CovMat.npy'%(Sims_Filename))
	errCFm_c_MUM_Corrected = temp.diagonal()*Rescale
	errCF_c_MUM_Corrected = np.vstack(( errCFp_c_MUM_Corrected, errCFm_c_MUM_Corrected ))

	temp = np.load('%s.SS112.rCLIP_X1sigma.MaskedCxUnmaskedC+.CovMat.npy'%(Sims_Filename))
	errCFp_c_MUM = temp.diagonal()*Rescale
	temp = np.load('%s.SS112.rCLIP_X1sigma.MaskedCxUnmaskedC-.CovMat.npy'%(Sims_Filename))
	errCFm_c_MUM = temp.diagonal()*Rescale
	errCF_c_MUM = np.vstack(( errCFp_c_MUM, errCFm_c_MUM ))



	# The pure noise
	Sims_Filename = '%s/Correlation_Function/100Sqdeg_NOISE_W3Mask_8.53GpAM_zKiDS_ZBcutNone/ThBins%s/NLOS114/NOISE_W3Mask.8.53GpAM.NLOS114' %(overall_DIR,ThBins_Sims)
	errCFp_c4 = np.loadtxt('%s.SS112.rCLIP_X1sigma.AverageCF+.asc' %Sims_Filename, usecols=(2,), unpack=True)
	errCFm_c4 = np.loadtxt('%s.SS112.rCLIP_X1sigma.AverageCF-.asc' %Sims_Filename, usecols=(2,), unpack=True)	
	errCF_c4 = np.vstack(( errCFp_c*np.sqrt(Rescale), errCFm_c*np.sqrt(Rescale) ))
	errCF_c5 = errCF_c4

	# The corrected errors
	errCFp_c_Corrected = np.loadtxt('%s.SS112.rCLIP_X1sigma.AverageCF+.Corrected.asc' %Sims_Filename, usecols=(2,), unpack=True)
	errCFm_c_Corrected = np.loadtxt('%s.SS112.rCLIP_X1sigma.AverageCF-.Corrected.asc' %Sims_Filename, usecols=(2,), unpack=True)	
	errCF_c_Corrected = np.vstack(( errCFp_c_Corrected*np.sqrt(Rescale), errCFm_c_Corrected*np.sqrt(Rescale) ))
	NLOS=1 # Redefine to avoid any more rescaling







# Save the corrected clipped CF and its error
np.savetxt('%s.SS%s.rCLIP_%ssigma.AverageCF+.Corrected.asc'%(combined_name,SS,sigma), np.c_[theta, CF_c_Corrected[0,:], errCF_c_Corrected[0,:]], header = 'theta/arcmin mean(xi_p) err(xi_p)', comments='#')

np.savetxt('%s.SS%s.rCLIP_%ssigma.AverageCF-.Corrected.asc'%(combined_name,SS,sigma), np.c_[theta, CF_c_Corrected[1,:], errCF_c_Corrected[1,:]], header = 'theta/arcmin mean(xi_m) err(xi_m)', comments='#')




def Calc_errRatio(num, errnum, denom, errdenom, crosscov):
	Ratio = num / denom
	errRatio = abs(Ratio)*np.sqrt( abs((errnum/num)**2. + (errdenom/denom)**2. -
 2.*(crosscov/(denom*num))) ) 
	return (Ratio - 1.), errRatio # return % difference



pm=['+', '-']


for i in range(len(pm)):

	plt.figure()
	plt.xscale('log')
	plt.errorbar(theta, theta*CF_uc[i,:]*1.e4, yerr=theta*errCF_uc[i,:]*1.e4/np.sqrt(float(NLOS)), color='dimgrey', linewidth=3.0, label = r'Unclipped, mask')
	plt.errorbar(np.exp(np.log(theta)+0.05), theta*CF_uc2[i,:]*1.e4, yerr=theta*errCF_uc2[i,:]*1.e4/np.sqrt(float(NLOS)), color='darkblue', linewidth=3.0, label = r'Unclipped, no mask')

	plt.errorbar(np.exp(np.log(theta)-0.01), theta*CF_c2[i,:]*1.e4, yerr=theta*errCF_c2[i,:]*1.e4/np.sqrt(float(NLOS)), color='cyan', linewidth=3.0, label = r'Clipped, no mask')
	#plt.errorbar(np.exp(np.log(theta)-0.05), theta*CF_c3[i,:]*1.e4, yerr=theta*errCF_c3[i,:]*1.e4/np.sqrt(float(NLOS)), color='lawngreen', linewidth=3.0, label = r'Clipped, W3 mask')	
	plt.errorbar(np.exp(np.log(theta)-0.05), theta*CF_c[i,:]*1.e4, yerr=theta*errCF_c[i,:]*1.e4/np.sqrt(float(NLOS)), color='magenta', linewidth=3.0, label = r'Clipped, mask')
	plt.errorbar(np.exp(np.log(theta)-0.05), theta*CF_c_Corrected[i,:]*1.e4, yerr=theta*errCF_c_Corrected[i,:]*1.e4/np.sqrt(float(NLOS)), color='purple', linewidth=3.0, label = r'Clipped, corrected')

	plt.errorbar(theta, theta*CF_c4[i,:]*1.e4, yerr=theta*errCF_c4[i,:]*1.e4/np.sqrt(float(NLOS)), color='orange', linewidth=3.0, label = r'Clipped, mask, NOISE')
	plt.errorbar(theta, theta*CF_c5[i,:]*1.e4, yerr=theta*errCF_c5[i,:]*1.e4/np.sqrt(float(NLOS)), color='red', linewidth=3.0, label = r'Clipped, no mask, NOISE')


	plt.xlim([0.8*np.min(theta), 1.2*np.max(theta)])
	#plt.ylim([0., 1.5*np.max(theta*CF_uc[i,:]*1.e4)])
	plt.xlabel(r'$\theta$ [arcmin]')
	plt.ylabel(r'$\theta \xi_{%s} \times 10^{-4}$'%pm[i])
	plt.legend(loc='best')
	plt.savefig('%s.SS%s.rCLIP_%ssigma.CF%s_MaskNoMaskComparison.eps' %(combined_name,SS,sigma,pm[i]))
	plt.show()



	# plot % diff between masked and unmasked. Assume for no cross-covariance.
	# G9 masked
	err=np.sqrt(float(NLOS))
	Ratio_uc, errRatio_uc = Calc_errRatio(CF_uc2[i,:], errCF_uc2[i,:]/err, CF_uc[i,:], errCF_uc[i,:]/err, errCF_uc_MUM[i,:]/err**2. )
	Ratio_c, errRatio_c = Calc_errRatio(CF_c2[i,:], errCF_c2[i,:]/err, CF_c[i,:], errCF_c[i,:]/err, errCF_c_MUM[i,:]/err**2.) 

	# W3 masked
	#Ratio_ucW3, errRatio_ucW3 = Calc_errRatio(CF_uc2[i,:], errCF_uc2[i,:]/err, CF_uc3[i,:], errCF_uc3[i,:]/err, errCF_uc_MUMW3/err**2. )
	#Ratio_cW3, errRatio_cW3 = Calc_errRatio(CF_c2[i,:], errCF_c2[i,:]/err, CF_c3[i,:], errCF_c3[i,:]/err, errCF_c_MUMW3/err**2.) 
	
	# Calculate error with no cross-covariance to compare
	Ratio_uc, errRatio_uc_NoCC = Calc_errRatio(CF_uc2[i,:], errCF_uc2[i,:]/err, CF_uc[i,:], errCF_uc[i,:]/err, 0.)
	Ratio_c, errRatio_c_NoCC = Calc_errRatio(CF_c2[i,:], errCF_c2[i,:]/err, CF_c[i,:], errCF_c[i,:]/err, 0.)



	# Correcting for the mask
	Ratio_c_Corrected, errRatio_c_Corrected = Calc_errRatio(CF_c2[i,:], errCF_c2[i,:]/err, CF_c_Corrected[i,:], errCF_c_Corrected[i,:]/err, errCF_c_MUM_Corrected[i,:]/err**2.)	

	

	plt.figure()
	plt.xscale('log')
	plt.errorbar(np.exp(np.log(theta)+0.05), Ratio_uc, yerr=errRatio_uc, color='dimgrey', linewidth=3.0, label = r'Unclipped, mask')
	#plt.errorbar(np.exp(np.log(theta)+0.01), Ratio_uc, yerr=errRatio_uc_NoCC, color='black', linewidth=3.0, label = r'Unclipped, No CC')

	#plt.errorbar(theta, Ratio_cW3, yerr=errRatio_cW3, color='lawngreen', linewidth=3.0, label = r'Clipped, W3 mask')
	plt.errorbar(np.exp(np.log(theta)-0.05), Ratio_c, yerr=errRatio_c, color='magenta', linewidth=3.0, label = r'Clipped, mask')
	plt.errorbar(theta, Ratio_c_Corrected, yerr=errRatio_c_Corrected, color='orange', linewidth=3.0, label = r'Clipped, mask Corrected')
		
	#plt.errorbar(np.exp(np.log(theta)+0.05), Ratio_c, yerr=errRatio_c_NoCC, color='orange', linewidth=3.0, label = r'Clipped, No CC')

	plt.plot( np.array([0.001, 500.]), np.array([-0.05, -0.05]), 'k--', linewidth=1.0)
	plt.plot( np.array([0.001, 500.]), np.array([0.05, 0.05]), 'k--', linewidth=1.0)

	plt.xlim([0.8*np.min(theta), 1.2*np.max(theta)])
	plt.ylim([ -0.2, 0.5])
	plt.xlabel(r'$\theta$ [arcmin]')
	plt.ylabel(r'$(\rm{NoMask} - \rm{Mask})/\rm{Mask}$ $\xi_{%s} $'%pm[i])
	plt.title(r'$\sigma_e \sim \mathcal{N}(0,0.29)$')
	plt.legend(loc='best')
	plt.savefig('%s.SS%s.rCLIP_%ssigma.CF%s_MaskNoMaskComparison_FracDiff.eps' %(combined_name,SS,sigma,pm[i]))
	plt.show()



	


 





