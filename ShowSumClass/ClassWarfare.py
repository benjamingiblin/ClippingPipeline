import sys
import numpy as np
import pylab as plt
from matplotlib import rc
rc('text',usetex=True)
rc('font',size=18)
rc('legend',**{'fontsize':18})
rc('font',**{'family':'serif','serif':['Computer Modern']})



class Filter_Input:
	def __init__(self, arguments):
		self.arguments = arguments

	def Filter(self):
		if self.arguments[1] == 'KiDS_Run':
			# Check if they're giving multiple fields as input (i.e. for plotting code)
			# Allow if so.
			Field=[]
			for f in self.arguments[3:]:
				if list(f)[0] == 'G':
					Field.append(f)
					#print("Identified target field %s" %f)
			try:
				los_end = int(self.arguments[-1])
			except ValueError:
				los_end=-1
			try:
				los = int(self.arguments[-2])
			except (ValueError, IndexError):
				los=-1 

			# Subtract off no. of args that are Field names.
			num_arg = 2 + len(Field)
			if los >= 0:
				num_arg += 2
				#print("Allowed no. of args is %s"%num_arg)
		

		elif self.arguments[1] == 'Sims_Run':
			num_arg = 4
		else:
			print("FIRST ARGUMENT MUST BE EITHER: \n")
			print("		KiDS_Run (for KiDS-450) \n")
			print("		Sims_Run (for 36sqdeg, 60sqdeg, 100sqdeg, OR 5000sqdeg(Mira Titan) mocks) \n")
			sys.exit(1)

		if len(self.arguments)-1 != num_arg:
			print('INVALID NUMBER OF ARGUMENTS.' )
			print('		KiDS arguments = KiDS_Run, paramfile, G.., G.., ..etc. Noise_los, Noise_los_end')
			print('							->Tag as many KiDS fields as desired on the end.')
			print('		Sims arguments = Sims_Run, paramfile, LOS_start, LOS_end')
			print('							->For Mira Titan run, LOS_start, LOS_end are dummies')
			print('PLEASE TRY AGAIN')
			sys.exit(1)	




	def Unpack_KiDS(self):
	
		# See if there are INTEGERS at the end of Field arguments
		# If so, these are the PURE NOISE REALISATION Numbers (used to correct for KiDS-450 mask)
		try:
			los_end = int(self.arguments[-1])
		except ValueError:
			los_end=-1
		try:
			los = int(self.arguments[-2])
		except (ValueError, IndexError):
			los=-1


		paramfile = self.arguments[2]
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

		Blind = args[0]
		SS=args[1]
		sigma=args[2]
		zlo=args[3]
		zhi=args[4]
		ThBins=args[5]
		MRres=args[6]
		OATH=args[7]
		NOISE=args[8]
		mask_variable=args[9]


		# Check if they're giving multiple fields as input (i.e. for plotting code)
		# Allow if so.
		Field=[]
		for f in self.arguments[3:]:
			if list(f)[0] == 'G':
				if NOISE=="Y" and los >= 0 and los_end >= 0:
					for nn in range(los, los_end+1): 
						Field.append('%s_NOISE%s'%(f,nn))
				else:
					Field.append(f)

		if len(Field) == 0:
			print("No valid KiDS fields given as input")
			sys.exit(1)
		
		

		# Check if there's a zB cut
		try:
			float(zlo)
			zlo_variable=True
		except ValueError:
			zlo_variable=False
			ZBcut='None'
		
		try: 
			float(zhi)
			zhi_variable=True		
		except ValueError:
			zhi_variable=False
			ZBcut='None'

		if zlo_variable and zhi_variable:
			if zhi > zlo:
				ZBcut='%s-%s' %(zlo, zhi)
			else:
				ZBcut='None'

		if MRres == '' or MRres == '-' or MRres == 'arcmin' or MRres == '1arcmin':
			MRres='arcmin'
			Prepend=''
		else:
			Prepend='MRres%s_' %MRres


		# Check if we're running shape noise only on this KiDS_Run
		if NOISE == "N" or NOISE == "" or NOISE == "-":
			NOISE = "N"
			Add_Noise = ""
		else:
			NOISE = "Y"
			Add_Noise = "_NOISE" 


		# Check if we're doing masking on this KiDS_Run 
		if mask_variable == "mask" or mask_variable == "" or mask_variable == "-":
			mask_variable=""
			Add_Mask=""
		else:
			mask_variable="-nomask"
			Add_Mask="_NoMask"


		DIRname = '%sKiDS_Fields_ZBcut%s%s%s' %(Prepend, ZBcut, Add_Noise, Add_Mask)

			
		# Check if No. Athena Theta bins is an integer
		try:
			int(ThBins)
		except ValueError:
			print("The number of Athena theta bins in param file is not an integer. \n Setting this number to default of 9 bins")
			ThBins=9
			
		return DIRname, Blind, SS, sigma, zlo, zhi, ThBins, OATH, Field[0]



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


		# Check if shape noise
		if SN == 'ALL' or SN == 'All' or SN == 'all':
			name_start = 'NOISE_'
			DIRname_start='_NOISE'
		elif 'Cycle' in SN or 'cycle' in SN:
			if float(gpam) == 3.32:
				name_start='SN0.28_'
			else:
				name_start='SN0.29_'
			DIRname_start='_SNCycle'
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
		else:
			print("mask variable in paramfile must be 'mask' or 'nomask'. Please fix this.")
			sys.exit(1)

		# Check if high/low sigma 8
		if cosmol == 'high_sigma8':
			DIRname_end='_HS8'
		elif cosmol == 'low_sigma8':
			DIRname_end='_LS8'
		elif cosmol == 'fid':
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
		elif int(sqdeg) == 5000 :
			Prepend='NSIDE%s_'%MRres 	# 'NSIDExxxx_'
		

		name = '%s%s' %(name_start, name_end)
		# e.g. NOISE_Mask, NF_test, SN0.27_Mask etc.
		DIRname = '%s%sSqdeg%s%s%sGpAM_z%s_ZBcut%s%s' %(Prepend,sqdeg, DIRname_start, DIRname_mid, gpam, z, ZBcut, DIRname_end) 
		z='%s%s' %(z, DIRname_end)
		print("DIRname is %s, z is %s" %(DIRname, z))
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
		plt.xscale('log')
		plt.yscale('log')
		plt.imshow(CCC_Mat, vmin=CCC_Mat.min(), vmax=CCC_Mat.max(), origin='lower', interpolation='nearest', extent = [theta.min(), theta.max(), theta.min(), theta.max()])
		plt.ylabel(r'$\theta$ [arcmin]')
		plt.xlabel(r'$\theta^{\prime}$ [arcmin]')
		plt.colorbar()
		plt.savefig('%sCCC_Mat.eps'%(savename))
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
	def Plot_CFs(self, theta_array, CF_array, errCF_array, legend_array, colour_array, theta_theory, CF_theory, savename, pm):

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
		theta_array = np.reshape(theta_array, (no_CFs, self.no_bins) )



		plt.figure()
		plt.xscale('log')
		plt.plot(theta_theory, theta_theory*CF_theory*1.e4, 'b:', linewidth=3.0, label = r'Takahashi$+12$')

		for i in range(0, len(theta_array[:,0])):

			eb1=plt.errorbar(theta_array[i,:], theta_array[i,:]*CF_array[i,:]*1.e4, yerr=theta_array[i,:]*errCF_array[i,:]*1.e4, color=colour_array[i], linewidth=3.0, label = r'%s' %legend_array[i])
			#eb1[-1][0].set_linestyle('--') 


		plt.xlim([0.8*np.min(theta_array[:,0]), 1.2*np.max(theta_array[:,-1])])
		plt.ylim([0., 1.5*np.max(theta_theory*CF_theory*1.e4)])
		plt.xlabel(r'$\theta$ [arcmin]')
		plt.ylabel(r'$\theta \xi_{%s} \times 10^{-4}$'%pm)
		plt.legend(loc='best')
		plt.savefig('%s.CF%s.eps'%(savename,pm))
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
		plt.savefig('%s.CF%sRatio.eps'%(savename,pm))
		#plt.show()	

		return




				
				

			
































		
		







