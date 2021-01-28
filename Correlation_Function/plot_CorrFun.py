# 26/10/2015
# B. M. Giblin, PhD student EDinburgh
# Calculate the average xi+/xi- CFs over NLOS, save and plot.
# AVERAGE over NLOS number of LOS


import pylab as plt
import numpy as np
import sys
import os
import matplotlib
from subprocess import call
import glob
from os import getcwd
from matplotlib import rc
from matplotlib import rcParams
from scipy.interpolate import interp1d

# Some font setting
rcParams['ps.useafm'] = True
rcParams['pdf.use14corefonts'] = True

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)
plt.rcParams["mathtext.fontset"] = "cm"

pipeline_DIR='/home/bengib/Clipping_SimsLvW/' # Where the codes of the pipeline reside
data_DIR='/home/bengib/Clipping_SimsLvW/' # Where the outputfiles and data are saved mid-computation
classdir = pipeline_DIR + "/ShowSumClass"

sys.path.insert(0, classdir) # add directory in which classes & functions 
							 # are defined to the python path
from ClassWarfare import Filter_Input, Handle_CF_Files
from FunkShins import Sort_Array_IntoGroups

variable = Filter_Input(sys.argv)
variable.Filter()
RUN = sys.argv[1]




TC='Y' 					# If TC='Y', then this will average/plot TreeCorr instead of Athena
Components_Calc = "N"	# If the delta-shear and deltaXORIG-shear files exist, AND it's a Sims_Run, this will be changed.


if TC == 'Y':
	Tree='Tree_'
	theta_col=1
	xip_col=3 # xip is in a different place in TreeCorr output file to Athena
	w_col=9
else:
	Tree=''
	theta_col=0
	xip_col=1
	w_col=7

if RUN == 'Sims_Run':

	name, gpam, DIRname, SS, sigma, SN, mask, z, PS, sqdeg, zlo, zhi, ThBins, OATH, los_start, los_end = variable.Unpack_Sims()
	
	indir='%s/%sCorrelation_Function/%s/ThBins%s/' %(data_DIR, Tree, DIRname, ThBins)
	# Group the input files
	files_uc = glob.glob('%s%s.%sGpAM.LOS*.ORIG.CorrFun.asc'%(indir,name,gpam))
	files_c = glob.glob('%s%s.%sGpAM.LOS*.SS%s.rCLIP_%ssigma.CorrFun.asc'%(indir,name,gpam,SS,sigma))
	frac_files = glob.glob('%s/Correlation_Function/%s/Clipped_PxlFrac.%s.%sGpAM.LOS*.SS%s.%ssigma.txt' %(pipeline_DIR, DIRname, name, gpam, SS, sigma))
	vol_files = glob.glob('%s/Correlation_Function/%s/Clipped_PxlVol.%s.%sGpAM.LOS*.SS%s.%ssigma.txt' %(pipeline_DIR, DIRname, name, gpam, SS, sigma))
	# Remove files outside of the LOS range.
	files_uc_copy = []
	files_c_copy = []
	frac_files_copy = []
	vol_files_copy = []

	for f in files_uc:
		try:
			l = int(f.split('LOS')[-1].split('.ORIG')[0])
		except ValueError:
			l = int(f.split('LOS')[-1].split('.ORIG')[0].split('n')[0])
		if l <= int(los_end) and l >= int(los_start):
			files_uc_copy.append(f)

	for i in range(len(files_c)):
		try:
			l = int(files_c[i].split('LOS')[-1].split('.SS')[0])
		except ValueError:
			l = int(files_c[i].split('LOS')[-1].split('.SS')[0].split('n')[0]) 
		if l <= int(los_end) and l >= int(los_start):
			files_c_copy.append(files_c[i])
			frac_files_copy.append(frac_files[i])	
			vol_files_copy.append(vol_files[i])
	files_uc = files_uc_copy
	files_c = files_c_copy
	frac_files = frac_files_copy
	vol_files = vol_files_copy

	# Average also the delta-shear and deltaXORIG-shear if these files exist
	if os.path.isfile('%s%s.%sGpAM.LOS%s.SS%s.deltaCLIP_%ssigma.CorrFun.asc'%(indir,name,gpam,los_start,SS,sigma)):
		Components_Calc = "Y"
		files_delta = glob.glob('%s%s.%sGpAM.LOS*.SS%s.deltaCLIP_%ssigma.CorrFun.asc'%(indir,name,gpam,SS,sigma))		
		files_deltaXORIG = glob.glob('%s%s.%sGpAM.LOS*.SS%s.deltaXORIGCLIP_%ssigma.CorrFun.asc'%(indir,name,gpam,SS,sigma))
		files_delta_copy = []
		files_deltaXORIG_copy = []
		for i in range(len(files_delta)):
			l = int(files_delta[i].split('LOS')[-1].split('.SS')[0])
			if l <= int(los_end) and l >= int(los_start):
				files_delta_copy.append(files_delta[i])
				files_deltaXORIG_copy.append(files_deltaXORIG[i])
		
		files_delta = files_delta_copy
		files_deltaXORIG = files_deltaXORIG_copy

	NLOS = len(files_uc)
	if os.path.isdir('%s/NLOS%s'%(indir, NLOS)) is False:
		call(['mkdir', '%s/NLOS%s'%(indir, NLOS)])
	outdir='%s/NLOS%s/'%(indir, NLOS)

	combined_name = '%s%s.%sGpAM.NLOS%.0f'%(outdir,name,gpam,NLOS)


	print('You have selected \n TreeCorr?: %s, \n Mock run: %s, \n Shape Noise: %s, \n Masking: %s, \n Gal. density: %s, \n F-Scale: %s Pxls, \n Clip threshold: %s, \n n(z): %s, \n zlow: %s, \n zhigh %s, \n Athena bins: %s, \n Output filename: %s, \n Output Subdir: %s, \n LOS start: %s, \n LOS end: %s, \n' %(TC, sqdeg, SN, mask, gpam, SS, sigma, z, zlo, zhi, ThBins, name, DIRname, los_start, los_end) )


else:
	
	DIRname, Blind, SS, sigma, zlo, zhi, ThBins, OATH, Dummy = variable.Unpack_KiDS()
	# check if multiple fields are given as input
	Field=[]
	arguments = sys.argv[3:]
	try:
		los_end = int(arguments[-1])
	except ValueError:
		los_end = -1
	try:
		los = int(arguments[-2])
	except (ValueError, IndexError):
		los = -1
	for f in arguments:
		if list(f)[0] == 'G':
			if "NOISE" in DIRname and los >= 0 and los_end >= 0:
				for nn in range(los, los_end+1): 
					Field.append('%s_NOISE%s'%(f,nn))
			else:
				Field.append(f)


	#arguments = sys.argv[3:]
	#Field = []
	#for f in arguments:
	#	if list(f)[0] == 'G':
	#		Field.append(f)

	indir='%s/%sCorrelation_Function/%s/ThBins%s/' %(data_DIR, Tree, DIRname, ThBins)
	# Group the input files and print some stuff to screen

	files_uc=[]
	files_c=[]
	frac_files=[]
	vol_files=[]
	filename='' # the beginning part of output filenames
				# ='G*G*G*...' depending on number of fields
	Field_Count=0 # count how many fields you ran on... use this to determine the no. of noise realisations
	print('You have selected Field(s):')
	for f in Field:
		print('%s '%f)
		files_uc.append('%s/%s.Blind%s.ORIG.CorrFun.asc' %(indir,f,Blind))
		files_c.append('%s/%s.Blind%s.SS%s.rCLIP_%ssigma.CorrFun.asc' %(indir,f,Blind,SS,sigma))
		frac_files.append('%s/Correlation_Function/%s/Clipped_PxlFrac.%s.Blind%s.SS%s.%ssigma.txt' %(pipeline_DIR,DIRname,f,Blind,SS,sigma))
		vol_files.append('%s/Correlation_Function/%s/Clipped_PxlVol.%s.Blind%s.SS%s.%ssigma.txt' %(pipeline_DIR,DIRname,f,Blind,SS,sigma))

		strip_name = f.split('_NOISE')[0]
		if strip_name not in filename:
			filename = filename + strip_name
			Field_Count+=1 # keep counting the fields thtat went int it
	if 'NOISE' in DIRname:
		filename = filename + '_NOISE%s' %int(len(Field)/Field_Count)

	print('TreeCorr? %s \nBlind %s \nF-scale %s Pxls, \nclip threshold %s sigma, \nzlow: %s, \nzhigh: %s, \nTheta Bins: %s' %(TC,Blind,SS,sigma,zlo,zhi,ThBins))
	

	outdir=indir

	combined_name = '%s/%s.Blind%s'%(outdir, filename, Blind)
	NLOS=len(Field)
	z='KiDSData'
		


no_bins = int(ThBins)
Clip_Class = Handle_CF_Files(files_c, no_bins)
Unclip_Class = Handle_CF_Files(files_uc, no_bins)
# Average the components as well if they exist.
if Components_Calc == "Y":
	Delta_Class = Handle_CF_Files(files_delta, no_bins)
	DeltaXUnclip_Class = Handle_CF_Files(files_deltaXORIG, no_bins)

	theta, CFp_d, theta_array, CFp_d_array, npair_array = Delta_Class.Average_CFs(theta_col, xip_col, w_col)
	theta, CFm_d, theta_array, CFm_d_array, npair_array = Delta_Class.Average_CFs(theta_col, xip_col+1, w_col)

	theta, CFp_dXO, theta_array, CFp_dXO_array, npair_array = DeltaXUnclip_Class.Average_CFs(theta_col, xip_col, w_col)
	theta, CFm_dXO, theta_array, CFm_dXO_array, npair_array = DeltaXUnclip_Class.Average_CFs(theta_col, xip_col+1, w_col)


# Average the CFs for each Field/los

theta, CFp_c, theta_array, CFp_c_array, npair_array = Clip_Class.Average_CFs(theta_col, xip_col, w_col)
theta, CFm_c, theta_array, CFm_c_array, npair_array = Clip_Class.Average_CFs(theta_col, xip_col+1, w_col)
									# column 1 = xi+, col 2 = xi-
theta, CFp_uc, theta_array, CFp_uc_array, npair_array = Unclip_Class.Average_CFs(theta_col, xip_col, w_col)
theta, CFm_uc, theta_array, CFm_uc_array, npair_array = Unclip_Class.Average_CFs(theta_col, xip_col+1, w_col)





# Find the bins and LOS that have anticorrelated clipped and unclipped
ij_array = np.empty([2])
antifiles = np.empty([])
for i in range(0, len(files_c)):
	for j in range(0, no_bins):
		if (CFp_uc_array[i,j] - CFp_uc[j]) * (CFp_c_array[i,j] - CFp_c[j]) < 0.:
			ij_array = np.vstack([ ij_array, np.array([i,j]) ])
			antifiles = np.append(antifiles, files_uc[i])
		
ij_array = np.delete(ij_array, 0, axis=0)




if NLOS > 1:

# CALCULATING COVARIANCE OF C/UC RATIOS IN 2 WAYS
# -----------------------------------------------------------
# 1. One in which you do the sum over each LOS.
#			This method is not great as we're taking the ratio of noisy numbers.
# 2. Instead, group the LOS into groups of 4 (each ~ area of KiDS-450),
#			Take the mean of each group of 4; these will be less noisy than the mean 
# 			ratio of each indiv LOS. Then take the covariance across these new means. 

	
	# Method 2.	
	Grouped_CFp_uc_array = Sort_Array_IntoGroups(CFp_uc_array, 4)
	Grouped_CFm_uc_array = Sort_Array_IntoGroups(CFm_uc_array, 4)
	Grouped_CFp_c_array = Sort_Array_IntoGroups(CFp_c_array, 4)
	Grouped_CFm_c_array = Sort_Array_IntoGroups(CFm_c_array, 4)
	
	# Mix of Methods 1. and 2. here
	Grouped_Ratio_CFp = np.zeros(no_bins)
	Grouped_Ratio_CFm = np.zeros(no_bins)
	Ratio_CFp = np.zeros(no_bins)
	Ratio_CFm = np.copy( Ratio_CFp )
	for i in range(0, no_bins):
		# Method 2.
		Grouped_Ratio_CFp[i] = np.mean( Grouped_CFp_c_array[:,i] / Grouped_CFp_uc_array[:,i] )
		Grouped_Ratio_CFm[i] = np.mean( Grouped_CFm_c_array[:,i] / Grouped_CFm_uc_array[:,i] )

		# Method 1.
		Ratio_CFp[i] = np.mean( CFp_c_array[:,i] / CFp_uc_array[:,i] )
		Ratio_CFm[i] = np.mean( CFm_c_array[:,i] / CFm_uc_array[:,i] )

	
	# Calculate the covariance --> Need to SQRT the returned error.

	# CLIPPED
	errCFp_c, Cov_CFp_c, CCC_CFp_c = Clip_Class.Calc_Covariance(CFp_c_array, CFp_c, CFp_c_array, CFp_c, '%s.SS%s.rCLIP_%ssigma.CxC+.'%(combined_name,SS,sigma))
	errCFm_c, Cov_CFm_c, CCC_CFm_c = Clip_Class.Calc_Covariance(CFm_c_array, CFm_c, CFm_c_array, CFm_c, '%s.SS%s.rCLIP_%ssigma.CxC-.'%(combined_name,SS,sigma))

	# UNCLIPPED
	errCFp_uc, Cov_CFp_uc, CCC_CFp_uc = Unclip_Class.Calc_Covariance(CFp_uc_array, CFp_uc, CFp_uc_array, CFp_uc, '%s.UCxUC+.'%(combined_name))
	errCFm_uc, Cov_CFm_uc, CCC_CFm_uc = Unclip_Class.Calc_Covariance(CFm_uc_array, CFm_uc, CFm_uc_array, CFm_uc, '%s.UCxUC-.'%(combined_name))

	# Unclipped X Clipped
	errCFp_cuc, Cov_CFp_cuc, CCC_CFp_cuc = Clip_Class.Calc_Covariance(CFp_c_array, CFp_c, CFp_uc_array, CFp_uc, '%s.SS%s.rCLIP_%ssigma.CxUC+.'%(combined_name,SS,sigma))
	errCFm_cuc, Cov_CFm_cuc, CCC_CFm_cuc = Clip_Class.Calc_Covariance(CFm_c_array, CFm_c, CFm_uc_array, CFm_uc, '%s.SS%s.rCLIP_%ssigma.CxUC-.'%(combined_name,SS,sigma))


	# Cov of Clip/Unclip ratio

	# Method 1.
	errCFp_Ratio, Cov_CFp_Ratio, CCC_CFp_Ratio = Clip_Class.Calc_Covariance(CFp_c_array/CFp_uc_array, Ratio_CFp, CFp_c_array/CFp_uc_array, Ratio_CFp, '%s.SS%s.rCLIP_%ssigma.Ratio+.'%(combined_name,SS,sigma))
	errCFm_Ratio, Cov_CFm_Ratio, CCC_CFm_Ratio = Clip_Class.Calc_Covariance(CFm_c_array/CFm_uc_array, Ratio_CFm, CFm_c_array/CFm_uc_array, Ratio_CFm, '%s.SS%s.rCLIP_%ssigma.Ratio-.'%(combined_name,SS,sigma))

	
	# Method 2. Breaks if you don't have enough CFs to sort into groups of 4.
	if NLOS > 8:
		errCFp_Ratio, Cov_CFp_Ratio, CCC_CFp_Ratio = Clip_Class.Calc_Covariance(Grouped_CFp_c_array/Grouped_CFp_uc_array, Grouped_Ratio_CFp, Grouped_CFp_c_array/Grouped_CFp_uc_array, Grouped_Ratio_CFp, '%s.SS%s.rCLIP_%ssigma.GroupedRatio+.'%(combined_name,SS,sigma))
		errCFm_Ratio, Cov_CFm_Ratio, CCC_CFm_Ratio = Clip_Class.Calc_Covariance(Grouped_CFm_c_array/Grouped_CFm_uc_array, Grouped_Ratio_CFm, Grouped_CFm_c_array/Grouped_CFm_uc_array, Grouped_Ratio_CFm, '%s.SS%s.rCLIP_%ssigma.GroupedRatio-.'%(combined_name,SS,sigma))



	# SQRT returned error
	errCFp_c = np.sqrt(errCFp_c)
	errCFm_c = np.sqrt(errCFm_c)

	errCFp_uc = np.sqrt(errCFp_uc)
	errCFm_uc = np.sqrt(errCFm_uc)

	# The cross-covariance term, DO NOT SQRT as it is sometimes -ve.



	# Plot the cross corr coeff matrices

	# CLIPPED
	#Clip_Class.Plot_Covariance(CCC_CFp_c, theta, '%s.SS%s.rCLIP_%ssigma.CxC+.'%(combined_name,SS,sigma))
	#Clip_Class.Plot_Covariance(CCC_CFm_c, theta, '%s.SS%s.rCLIP_%ssigma.CxC-.'%(combined_name,SS,sigma))
	# UNCLIPPED
	#Unclip_Class.Plot_Covariance(CCC_CFp_uc, theta, '%s.UCxUC+.'%(combined_name))
	#Unclip_Class.Plot_Covariance(CCC_CFm_uc, theta, '%s.UCxUC-.'%(combined_name))
	#UNCLIP x CLIP
	#Clip_Class.Plot_Covariance(CCC_CFp_cuc, theta, '%s.SS%s.rCLIP_%ssigma.CxUC+.'%(combined_name,SS,sigma))
	#Clip_Class.Plot_Covariance(CCC_CFm_cuc, theta, '%s.SS%s.rCLIP_%ssigma.CxUC-.'%(combined_name,SS,sigma))

else:
	
	errCFp_c = np.zeros(no_bins)
	errCFm_c = np.zeros(no_bins)
	errCFp_uc = np.zeros(no_bins)
	errCFm_uc = np.zeros(no_bins)

	errCFp_cuc = np.zeros(no_bins)
	errCFm_cuc = np.zeros(no_bins)




ZBcut = DIRname.split('ZBcut')[-1].split('_')[0]
# Save the average CFs: clipped/unclipped xi+/xi-
np.savetxt('%s.SS%s.rCLIP_%ssigma.AverageCF+.asc'%(combined_name,SS,sigma), np.c_[theta, CFp_c, errCFp_c], header = 'theta/arcmin mean(xi_p) std(xi_p)', comments='#')
np.savetxt('%s.SS%s.rCLIP_%ssigma.AverageCF-.asc'%(combined_name,SS,sigma), np.c_[theta, CFm_c, errCFm_c], header = 'theta/arcmin mean(xi_m) std(xi_m)', comments='#')

print('saving %s.AverageCF+.asc'%(combined_name))
np.savetxt('%s.AverageCF+.asc'%(combined_name), np.c_[theta, CFp_uc, errCFp_uc], header = 'theta/arcmin mean(xi_p) std(xi_p)', comments='#')
np.savetxt('%s.AverageCF-.asc'%(combined_name), np.c_[theta, CFm_uc, errCFm_uc], header = 'theta/arcmin mean(xi_m) std(xi_m)', comments='#')

if Components_Calc == "Y":
	np.savetxt('%s.SS%s.deltaCLIP_%ssigma.AverageCF+.asc'%(combined_name,SS,sigma), np.c_[theta, CFp_d, np.zeros(no_bins)], header = 'theta/arcmin mean(xi_p) std(xi_p)', comments='#')
	np.savetxt('%s.SS%s.deltaCLIP_%ssigma.AverageCF-.asc'%(combined_name,SS,sigma), np.c_[theta, CFm_d, np.zeros(no_bins)], header = 'theta/arcmin mean(xi_m) std(xi_p)', comments='#')
	
	np.savetxt('%s.SS%s.deltaXORIGCLIP_%ssigma.AverageCF+.asc'%(combined_name,SS,sigma), np.c_[theta, CFp_dXO, np.zeros(no_bins)], header = 'theta/arcmin mean(xi_p) std(xi_p)', comments='#')
	np.savetxt('%s.SS%s.deltaXORIGCLIP_%ssigma.AverageCF-.asc'%(combined_name,SS,sigma), np.c_[theta, CFm_dXO, np.zeros(no_bins)], header = 'theta/arcmin mean(xi_m) std(xi_p)', comments='#')



# Plot the CFs (stack the mean clipped & unclipped into array)
stack_theta = np.stack((theta, theta))
stack_CFp = np.stack((CFp_c, CFp_uc))
stack_CFm = np.stack((CFm_c, CFm_uc))

# similarly stack the errors on clipped and unclipped
stack_errCFp = np.stack((errCFp_c/np.sqrt(float(NLOS)), errCFp_uc/np.sqrt(float(NLOS)) ))
stack_errCFm = np.stack((errCFm_c/np.sqrt(float(NLOS)), errCFm_uc/np.sqrt(float(NLOS)) ))



legend_array = ['Clipped', 'Unclipped']
colour_array = ['magenta', 'darkblue']


# The plotting bit
#Clip_Class.Plot_CFs(stack_theta, stack_CFp, stack_errCFp, legend_array, colour_array, '%s.SS%s.rCLIP_%ssigma'%(combined_name,SS,sigma), '+')
#Clip_Class.Plot_CFs(stack_theta, stack_CFm, stack_errCFm, legend_array, colour_array, '%s.SS%s.rCLIP_%ssigma'%(combined_name,SS,sigma), '-')

stack_errCFp = np.stack((errCFp_c, np.zeros(len(theta)) ))
stack_errCFp_cuc = np.stack((errCFp_cuc, np.zeros(len(theta)) ))

stack_errCFm = np.stack((errCFm_c, np.zeros(len(theta)) ))
stack_errCFm_cuc = np.stack((errCFm_cuc, np.zeros(len(theta)) ))


#Clip_Class.Plot_CFs_Ratio(theta, stack_CFp, stack_errCFp, CFp_uc, errCFp_uc, stack_errCFp_cuc, legend_array, colour_array, '%s.SS%s.rCLIP_%ssigma'%(combined_name,SS,sigma), '+')

#Clip_Class.Plot_CFs_Ratio(theta, stack_CFm, stack_errCFm, CFm_uc, errCFm_uc, stack_errCFm_cuc, legend_array, colour_array, '%s.SS%s.rCLIP_%ssigma'%(combined_name,SS,sigma), '-')


# Calc the averaged clipped pxl frac and save the result.
Clip_Class.Calc_Clip_Frac(frac_files, '%s.SS%s.rCLIP_%ssigma.AvClipFrac.txt'%(combined_name,SS,sigma))
#Clip_Class.Calc_Clip_Frac(vol_files, '%s.SS%s.rCLIP_%ssigma.AvClipVol.txt'%(combined_name,SS,sigma))












































