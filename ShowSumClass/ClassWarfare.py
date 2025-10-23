import sys
import numpy as np
import pylab as plt
from matplotlib import rc
from matplotlib import rcParams
# Some font setting
rcParams['ps.useafm'] = True
rcParams['pdf.use14corefonts'] = True

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)
plt.rcParams["mathtext.fontset"] = "cm"



class Filter_Input:
	def __init__(self, arguments):
		self.arguments = arguments

	def Filter(self):

		if self.arguments[1] == 'Sims_Run' or self.arguments[1] == 'KiDS_Run':
			num_arg = 4
		else:
			print("FIRST ARGUMENT MUST BE EITHER: \n")
			print("		KiDS_Run (for KiDS-1000) \n")
			print("		Sims_Run (for 36sqdeg, 60sqdeg, 100sqdeg, OR 5000sqdeg(Mira Titan) mocks) \n")
			sys.exit(1)

		if len(self.arguments)-1 != num_arg:
			print('INVALID NUMBER OF ARGUMENTS.' )
			print('		arguments = Sims_Run/KiDS_Run, paramfile, LOS_start, LOS_end')
			print('							->For Mira Titan run, LOS_start, LOS_end are dummies')
			print('PLEASE TRY AGAIN')
			sys.exit(1)	



	def Unpack_Sims(self):
		paramfile = self.arguments[2]
		los = self.arguments[3]
		los_end = self.arguments[-1] # may be the same as los, if los_end is not
                
		f = open(paramfile)
		params = f.readlines()
		args=[]
		for line in params:
			try:
				foo = line.strip('\n').split()[0]
				if foo !='#' and foo != '':
					args.append(foo) 
			except (ValueError,IndexError):
				break
				# ignore blank lines at the end


		# Build general filename prefix and output directory from inputs:
		gpam = args[0]
		SS=args[1]
		sigma=args[2]
		SN=args[3]
		mask=args[4]
		z=args[5]
		PS=args[6]
		sqdeg=args[7]
		cosmol=args[8]
		zlo=args[9]
		zhi=args[10]
		ThBins=args[11]
		MRres=args[12]
		OATH=args[13]
		try:
			Sys=args[14] # Opt sys param (IA1.0, BaryON, etc.)
		except IndexError:
			Sys=""
                
		# Check if shape noise
		if 'ALL' in SN or 'All' in SN or 'all' in SN:
			if 'KiDS' in SN:
				name_start = 'NOISE-KiDS_'
				DIRname_start='_NOISE-KiDS'
			else:
				name_start = 'NOISE_'
				DIRname_start='_NOISE'
		elif 'Cycle' in SN or 'cycle' in SN:
			# Set the noise level (not specified in param_file if doing Cycle):
			if 'KiDS1000' in str(gpam):
				if float(zlo)==0.1 and float(zhi)==1.2:
					SN_level=0.265
				else:
					# Must set SN to values measured in each KiDS1000 bin.
					bin_edges = np.array([0.1, 0.3, 0.5, 0.7, 0.9, 1.2])
					sigma_e_values = [0.270, 0.258, 0.273, 0.254, 0.270]
					# Note: this only works for the 5 zbins used for cosmic shear.
					idx_sig = np.where(float(zlo) == bin_edges)[0][0]
					SN_level = sigma_e_values[idx_sig]
				name_start='SN%s_' %SN_level
                                
			else:
			        name_start='SN0.28_' 
                                        
			DIRname_start='_SNCycle'
		elif 'KiDS' in SN:
		        name_start='SNKiDS_'
		        DIRname_start='_SNKiDS'
		else:
			try:
				float(SN)
			except ValueError:
				print("The shape noise (SN) variable in the paramfile must be either 'ALL' or a numerical value. Please fix this.")
				sys.exit(1)

			if float(SN) == 0.:
				name_start='NF_'
				DIRname_start='_NF'
			else:
				name_start='SN%s_'%SN
				DIRname_start='_SN%s'%SN

		# Check if masking
		if mask == 'nomask':
			name_end='test'
			DIRname_mid='_NoMask_'
		elif mask == 'mask':
			name_end='Mask'
			DIRname_mid='_Mask_'
		elif mask == 'G9mask':
			name_end='G9Mask'
			DIRname_mid='_G9Mask_'
		elif mask == 'W3mask':
			name_end='W3Mask'
			DIRname_mid='_W3Mask_'
		elif 'mosaic' in mask or 'Mosaic' in mask:
			name_end='Mosaic'
			DIRname_mid='_Mosaic_'
		else:
			print("mask variable in paramfile is set to %s; not a valid option! Please fix this." %mask)
			sys.exit(1)

		# Check if high/low sigma 8
		if cosmol == 'high_sigma8':
			DIRname_end='_HS8'
		elif cosmol == 'low_sigma8':
			DIRname_end='_LS8'
		elif cosmol == 'fid' or cosmol=='KiDS1000':
			DIRname_end='_Cosmol%s' %cosmol
		else:
			try:
				int(cosmol)
				DIRname_end='_Cosmol%s' %cosmol
			except ValueError:
				DIRname_end=''

#		elif isinstance(cosmol,int) == True:
#			# If cosmol is an integer, then you're running DH10 mocks
#			DIRname_end='_Cosmol%s' %cosmol
#		else:
#			print("DIRname_end is...")
#			DIRname_end=''


		# Check if there's a zB cut
		try:
			float(zlo)
			zlo_variable=True
		except ValueError:
			zlo='None'
			zlo_variable=False
			ZBcut='None'
		
		try: 
			float(zhi)
			zhi_variable=True		
		except ValueError:
			zhi='None'
			zhi_variable=False
			ZBcut='None'

		if zlo_variable and zhi_variable:
			if zhi > zlo:
				ZBcut='%s-%s' %(zlo, zhi)
			else:
				ZBcut='None'



		# Check if No. Athena Theta bins is an integer
		try:
			int(ThBins)
		except ValueError:
			print("The number of Athena theta bins in param file is not an integer. \n Setting this number to default of 9 bins")
			ThBins=9


		if int(sqdeg) == 100 or int(sqdeg) == 60 or int(sqdeg) == 36:		
			if MRres == '5arcs' or MRres == '' or MRres == '-':
				Prepend=''
			else:
				Prepend='MRres%s_' %MRres
	
			if "IA" in Sys or "Bary" in Sys or "dz" in Sys or "SLC" in Sys:
				Prepend += '%s_' %Sys

		elif int(sqdeg) == 5000 :
			Prepend='NSIDE%s_'%MRres 	# 'NSIDExxxx_'
		

		name = '%s%s' %(name_start, name_end)
		# e.g. NOISE_Mask, NF_test, SN0.27_Mask etc.
		DIRname = '%s%sSqdeg%s%s%sGpAM_z%s_ZBcut%s%s' %(Prepend,sqdeg, DIRname_start, DIRname_mid, gpam, z, ZBcut, DIRname_end) 
		####if self.arguments[1] == 'KiDS_Run': DIRname = 'KiDS1000Data_%s' %DIRname; fi
		z='%s%s' %(z, DIRname_end)
		#print("DIRname is %s, z is %s" %(DIRname, z))
		#e.g. 	60Sqdeg_NF_Mask_8.53GpAM_z0.640
		#e.g. 	100Sqdeg_SN0.27_NoMask_8.53GpAM_zKiDS_HS8

		
		return name, gpam, DIRname, SS, sigma, SN, mask, z, PS, int(sqdeg), str(zlo), str(zhi), ThBins, OATH, los, los_end 

		
		




class Format_2Darray:
	
	def __init__(self, mapp):
		self.mapp = mapp

	# Pad_2Darray embeds a 2D array in a frame of given dimensions composed of 
	# a single number (again which you specify).
	def Pad_2Darray(self, new_sizeX, new_sizeY, number):
		noex_cols = int((new_sizeX - len(self.mapp[0,:]))/2.) # no. rows to add on each side
		noex_rows = int((new_sizeY - len(self.mapp[:,0]))/2.)

		extra_rows = np.zeros([noex_rows, len(self.mapp[0,:]) ]) + number
		self.mapp = np.r_[extra_rows, self.mapp] # add on the start
		self.mapp = np.r_[self.mapp, extra_rows] # add them on the end

		extra_cols = np.zeros([ len(self.mapp[:,0]), noex_cols]) + number
		self.mapp = np.c_[extra_cols, self.mapp] # add on the start
		self.mapp = np.c_[self.mapp, extra_cols] # add on the end

		# rounding of num extra row/cols to intg can mean end map is 1pxl too small:
		if self.mapp.shape[0] == new_sizeY-1:
			# need an extra row, add to the bottom
			self.mapp = np.r_[ self.mapp, np.zeros([1,self.mapp.shape[1]])+number ]
		if self.mapp.shape[1] == new_sizeX-1:
			# need an extra column, add it on the right side
			self.mapp = np.c_[ self.mapp, np.zeros([self.mapp.shape[0],1])+number ]
		        
		return self.mapp



	def Manipulate_W3(self, new_sizeX, new_sizeY, number):
		print("Manpulating the W3 mask")
		if new_sizeX < len(self.mapp[0,:]) and new_sizeY < len(self.mapp[:,0]):
			self.mapp = self.mapp[:new_sizeY, :new_sizeX]
		else:
			self.mapp = self.Pad_2Darray(new_sizeX, new_sizeY, 1.)
		return self.mapp


	# Do the opposite of Pad_2Darray --> cut map down to certain size
	def Cut_Edges(self, new_sizeX, new_sizeY):
	
		low_limX = int((len(self.mapp[0,:]) - new_sizeX)/2.)
		up_limX = low_limX + new_sizeX 

		low_limY = int((len(self.mapp[:,0]) - new_sizeY)/2.)
		up_limY = low_limY + new_sizeY 
		return self.mapp[low_limY:up_limY, low_limX:up_limX]













class Handle_CF_Files:
	
	def __init__(self, filelist, no_bins):
		self.filelist = filelist
		self.no_bins = no_bins

	
	def Unpack_CFs(self, column0,column,column2):
		# Unpack the files and return theta and CF function arrays.
		CF_array = np.zeros([1, self.no_bins])
		theta_array = np.zeros([1, self.no_bins])
		npair_array = np.zeros([1, self.no_bins])
		for f in self.filelist:
			theta, CF, npair = np.loadtxt(f, usecols=(column0,column,column2), unpack=True)
			theta_array = np.vstack([theta_array, theta])
			CF_array = np.vstack([CF_array, CF])
			npair_array = np.vstack([npair_array, npair])
		
		# remove first row of zeros
		theta_array = np.delete(theta_array, (0), axis=0)
		CF_array = np.delete(CF_array, (0), axis=0)
		npair_array = np.delete(npair_array, (0), axis=0)
		return theta_array, CF_array, npair_array




	def Average_CFs(self, column0,column,column2):
		# Take the average in each bin of the CF arrays
		theta_array, CF_array, npair_array = self.Unpack_CFs(column0,column,column2)
		meantheta = np.zeros(self.no_bins)
		meanCF = np.zeros(self.no_bins)
				
		for i in range(0, self.no_bins):
			meantheta[i] = np.mean(theta_array[:,i])
			if int(np.sum(npair_array[:,i])) != 0:
				meanCF[i] = np.sum(CF_array[:,i]*npair_array[:,i]) / np.sum(npair_array[:,i])
			else:
				meanCF[i] = 0.

		return meantheta, meanCF, theta_array, CF_array, npair_array





	def Calc_Covariance(self, CF_array1, meanCF1, CF_array2, meanCF2, savename):
		no_CFs = len(CF_array1[:,0])
		# cov matrix
		Cov_Mat = np.empty([self.no_bins, self.no_bins])
		for i in range(0, self.no_bins):
			for j in range(0, self.no_bins):
					Cov_Mat[i,j] = (1./(no_CFs-1.)) * np.sum( (CF_array1[:,i] - meanCF1[i]) * (CF_array2[:,j] - meanCF2[j]) )

		err = Cov_Mat.diagonal()
		# Save the cov matrix
		np.save('%sCovMat' %savename, Cov_Mat)
		
		# cross-corr coeff matrix
		CCC_Mat =  np.empty([self.no_bins, self.no_bins])
		for i in range(0, self.no_bins):
			for j in range(0, self.no_bins):		
				CCC_Mat[i,j] = Cov_Mat[i,j] / np.sqrt(  abs(Cov_Mat[i,i]*Cov_Mat[j,j]) )

		return err, Cov_Mat, CCC_Mat




	def Plot_Covariance(self, CCC_Mat, theta, savename):

		plt.figure()
		
		plt.imshow(CCC_Mat, vmin=CCC_Mat.min(), vmax=CCC_Mat.max(), origin='lower', interpolation='nearest', extent = [theta.min(), theta.max(), theta.min(), theta.max()])
		plt.ylabel(r'$\theta$ [arcmin]')
		plt.xlabel(r'$\theta^{\prime}$ [arcmin]')
		plt.colorbar()
		plt.savefig('%sCCC_Mat.png'%(savename))
		#plt.show()
		return



	def Calc_Clip_Frac(self, frac_files, savename):
		# Calculate the average clipped pxl fraction given a list of clipped pxl files
		# and the name under which to save the result.
		frac=[]
		for f in frac_files:
			thefile = open(f)
			content = thefile.readlines()
			frac.append(content[1].strip())

		meanfrac = np.mean(np.array(frac).astype(float))
		stdfrac = np.std(np.array(frac).astype(float))


		f = open(savename, 'w')
		f.write('The average fraction of pxls clipped across the kappa maps = \n')
		f.write('%s pm %s \n' %(meanfrac, stdfrac))
		f.close()

		print('The average fraction of pxls clipped across the kappa maps = %s pm %s' %(meanfrac, stdfrac))







 

	# last argument pm+ '+' or '-'
	def Plot_CFs(self, theta_array, CF_array, errCF_array, legend_array, colour_array, savename, pm):

		if pm != '+' and pm != '-':
			print('	ARGUMENTS:')
			print('		stacked theta array [arcmin]')
			print('		stacked CF array')
			print('		stacked error ON MEAN CF array')
			print('		stacked legend array')
			print('		stacked color_array')
			print('----------------------------')
			print('		theta theory [arcmin]')
			print('		CF theory [arcmin]')
			print('		savename')
			print('		+ OR -')
			sys.exit(1)
		 

		# Check if CF_array is 1D or 2D, and therefore get number of CFs to plot
		if len(CF_array.shape) ==1:
			# then it's 1D
			no_CFs = 1
		else:
			# then it's 2D
			no_CFs = len(CF_array[:,0])


		# Reshape the arrays: means code won't break for 1D CF_array
		CF_array = np.reshape(CF_array, (no_CFs, self.no_bins) )
		errCF_array = np.reshape(errCF_array, (no_CFs, self.no_bins) )
		

		plt.figure()
		plt.xscale('log')

		for i in range(0, len(theta_array[:,0])):

			eb1=plt.errorbar(theta_array[i,:], theta_array[i,:]*CF_array[i,:]*1.e4, yerr=theta_array[i,:]*errCF_array[i,:]*1.e4, color=colour_array[i], linewidth=3.0, label = r'%s' %legend_array[i])
			#eb1[-1][0].set_linestyle('--') 


		plt.xlim([0.8*np.min(theta_array[:,0]), 1.2*np.max(theta_array[:,-1])])
		plt.ylim([ np.min(theta_array[0,:]*CF_array[0,:]*1.e4),
                           1.5*np.max(theta_array[0,:]*CF_array[0,:]*1.e4) ])
		plt.xlabel(r'$\theta$ [arcmin]')
		plt.ylabel(r'$\theta \xi_{%s} \times 10^{-4}$'%pm)
		plt.legend(loc='best')
		plt.savefig('%s.CF%s.png'%(savename,pm))
		#plt.show()

		return







	def Plot_CFs_Ratio(self, theta, CF_array, errCF_array, CF_denom, errCF_denom, errCF_cc_array, legend_array, colour_array, savename,pm):


		if pm != '+' and pm != '-':
			print('	ARGUMENTS:')
			print('		theta [arcmin]')
			print('		stacked CFs you want to plot')
			print('		stacked error on CFs you want to plot [NOT ERROR ON MEAN]')
			print('		The CF you want on denominator of ratio')
			print('		The error on the CF you want on denominator of ratio')
			print('		stacked cross covariance between numerator CFs and denominator CF')
			print('		stacked legend array')
			print('		stacked color_array')
			print('----------------------------')
			print('		savename')
			print('		+ OR -')
			sys.exit(1)


		NLOS = len(self.filelist)

		# Check if CF_array is 1D or 2D, and therefore get number of CFs to plot
		if len(CF_array.shape) ==1:
			# then it's 1D
			no_CFs = 1
		else:
			# then it's 2D
			no_CFs = len(CF_array[:,0])


		# Reshape the arrays: means code won't break for 1D CF_array
		CF_array = np.reshape(CF_array, (no_CFs, self.no_bins) )
		errCF_array = np.reshape(errCF_array, (no_CFs, self.no_bins) )
		errCF_cc_array = np.reshape(errCF_cc_array, (no_CFs, self.no_bins) )


		Ratios = np.zeros([no_CFs, self.no_bins ])
		errRatios = np.zeros([no_CFs, self.no_bins ])
		for i in range(0, no_CFs):
			Ratios[i,:] = CF_array[i,:]/CF_denom
			errRatios[i,:] = abs(Ratios[i,:])*np.sqrt( (errCF_array[i,:]/CF_array[i,:])**2. + (errCF_denom/CF_denom)**2. - 2.*(errCF_cc_array[i,:]/(CF_denom*CF_array[i,:])) ) / NLOS

			Ratios[i, np.where( np.isfinite(Ratios[i,:]) == False ) ] = 0.
			errRatios[i, np.where( np.isfinite(errRatios[i,:]) == False ) ] = 0.

		# Do not square the cross-covariance term. It is already squared. 
		# (i.e. I never took the sqrt as it is sometimes -ve). 		
		# Division by self.filelist makes it error on MEAN


		# Plot the ratio 
		plt.figure()
		plt.xscale('log')
		for i in range(0, no_CFs):
			eb1=plt.errorbar(theta, Ratios[i,:], yerr=errRatios[i,:], color=colour_array[i], linewidth=3.0, label = r'%s' %legend_array[i])


		plt.xlim([0.8*np.min(theta[0]), 1.2*np.max(theta[-1])])
		#plt.ylim([0.8*np.min(Ratios), 1.3*np.max(Ratios)])
		plt.ylim([0.2, 1.8])
		plt.xlabel(r'$\theta$ [arcmin]')
		plt.ylabel(r'$\xi_{%s}^{c}/\xi_{%s}^{uc}$'%(pm,pm))
		plt.legend(loc='best')
		plt.savefig('%s.CF%sRatio.png'%(savename,pm))
		#plt.show()	

		return




				
				

			
































		
		







