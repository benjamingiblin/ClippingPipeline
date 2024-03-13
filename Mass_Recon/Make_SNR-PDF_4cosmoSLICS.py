# 06/11/2019, B. M. Giblin, Postdoc, Edinburgh
# Read in the Kappa Maps for the cosmoSLICS cosmologies, convert to PDFs and save.

import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import sys
import time
from scipy.interpolate import interp1d
from astropy.io import fits
from natsort import natsorted # used to sort strings in numerical order
                              # (i.e. it knows '2' comes before '10') 

def View_PDF(bin_cen, pdf, pdf_avg, pdf_err): 
             #bin_cen2, pdf2, pdf_avg2, pdf_err2):
        plt.figure()
        for i in range(pdf.shape[0]):
                plt.plot(bin_cen, pdf[i,:], color='dimgrey')
        #for i in range(pdf2.shape[0]):
        #        plt.plot(bin_cen2, pdf2[i,:], color='yellow')
        plt.errorbar(bin_cen, pdf_avg, yerr=pdf_err, color='black', linewidth=3, label='orig')
        #plt.errorbar(bin_cen2, pdf_avg2, yerr=pdf_err2, color='orange', linewidth=3, label='re-binned')
        #plt.yscale('log')
        plt.xlabel('SNR')
        plt.ylabel('PDF(SNR)')
        plt.legend(loc='best', frameon=False)
        plt.show()
        return

def Remove_Mask_Spike_From_PDF(kappa, pdf):
        # Even when we only include finite pxls (avoiding masks)
        # pxls close to masked regions get biased low, making a spike
        # in the centre of the PDF. This removes it.
        i1 = int((nbins/2)-2) # bracket the central peak
        i2 = int((nbins/2)+2) # 2 idx above, 2 below
        tmp_pdf = np.append( pdf[:i1], pdf[i2:] ) # rm middle
        tmp_kap = np.append( kappa[:i1], kappa[i2:] )
        tmp_interp = interp1d( tmp_kap, tmp_pdf, kind='cubic') # interp across the gap
        new_pdf = tmp_interp( kappa ) # same as old pdf for all idx except central 4 idx.
        
        # old redundant way.
        #pdf[ np.argmax(pdf) ] = pdf[ np.argmax(pdf)-1 ]
        # (set peak value, at zero due to the mask) to value next to the peak)
        return new_pdf

def Read_AndOr_Make_Overall_Noise_Map(noise_fname, filenames):
        # check if the STDEV of all the noise maps is saved to file.                                                 
        # If not, read in all the noise maps and make the stdev map.
        if os.path.isfile(noise_fname):
                print("Reading in a pre-saved NOISE-STDEV map.")
                Noise_std = np.load(noise_fname)
        else:
                t1 = time.time()
                npix_x = np.load(filenames[0]).shape[1]
                npix_y = np.load(filenames[0]).shape[0]
                Noise_maps = np.zeros([ len(filenames), npix_y, npix_x ]) # Store all of the Noise maps 
                for i in range(len(filenames)):
                        #print( "Reading in NOISE file %s of %s"%(i,len(filenames_n)) )
                        Noise_maps[i] = np.load(filenames[i])
                        #small_pxls = len(np.where( np.logical_and(Noise_maps[i]>-0.001,
                                                                  #Noise_maps[i]<0.001) )[0])
                        # NB: previously used +/-0.005, but too wide for LSST-18.66-noise.
                        #small_pxls_pc = 100.*small_pxls / (npix_x * npix_y)
                        #if small_pxls_pc > 95:
                        #        print("!!! The following NOISE map has %.3f pecent small pxls !!!" %(small_pxls_pc))
                        #        print( filenames[i] )
                        #        print(" +++ DELETING IT +++ ")
                                #os.remove( filenames_n[i] )
                        #        print( " ----------------------------------------- " )
                                
                Noise_std = np.std( Noise_maps, axis=0 ) # 2D stdev map                             
                np.save( noise_fname, Noise_std ) #np.array([Noise_std]) )
                t2 = time.time()
                print("It took %.0f seconds to read in %s Noise maps and save the stdev map." %(t2-t1, len(filenames)) )
        return Noise_std


def Rebin_PDF(fine_bin_centres, fine_PDF, coarse_bin_edges):
        # Re-bin a finely-binned PDF into  coarsely-binned PDF
        coarse_PDF = np.zeros( nbins_coarse )
        idx = []
        for j in range(nbins_coarse):
                # find all the fine bins within this coarse bin
                id = np.where( np.logical_and( fine_bin_centres<(coarse_bin_edges[j+1]),
                                               fine_bin_centres>=(coarse_bin_edges[j]) ) )[0]
                idx+=[list(id)]
                coarse_PDF[j] = np.sum( fine_PDF[idx[j]] ) / len(fine_bin_centres)
        return coarse_PDF, idx


def Shuffle_filenames_LOSnR( filenames, R ):
        # Re-shuffle the filenames in accordance with the 217x18 shuffle matrix                      
        files_shffld = []
        for i in range(len(filenames_s)):
                los_num = filenames[i].split('LOS')[-1].split('R')[0]
                # replace the LOS num with another from the matrix
                # (which is different for each of the 18 regions):
                files_shffld.append( filenames[i].replace('LOS%sR%s' %(los_num,R),
                                                          'LOS%sR%s'%(int(shffld_lr[i,R-1]),R) ))
        return files_shffld


# ------------------------------------------------------------------------ #

# Define binning for the PDF
# If "Fine" it reads in a map and makes a finely-binned PDF
# If "Coarse" it reads in the finely-binned PDF and re-bins it,
# to make a coarsely-binned PDF.
# This saves us storing all the kappa maps (several TB).                                                                                
Binning_Fine_Coarse = "Coarse"

# ------------------------------------------------------------------------ #   

mock_Type = 'cosmoSLICS'    # 'cosmoSLICS' 'SLICS'
Survey = 'KiDS1000'    # 'LSST' 'KiDS1000'
Mask = 'Mosaic'        # "Mosaic' 'NoMask'
Read_IA = False        # if True read the measurements for the (fid-cosmol) IA-contam. measurements
A_IA = "1.0"           # IA amplitude

Calc_Combine_zbins = True   # If True, will read in the PDFs for combinations of redshift bins

Make_Data_Vec = False   # If True, 
                       # it will just read & avg the 1st 10 maps & avg them
                       # to make something with the noise levels of KiDS1000
Data_Vec_batch = 4     # The batch of LOS to avg (0=0:10, 1=10:20, 2=20:30, ...)
MRres = '140.64arcs'
SS = [2.816] # (1.408, 2.816, 5.631), (3.11, 9.33, 18.66, 84.85)

# binning scheme for the PDF
nbins_coarse = 4
nbins_fine = 2000

if Binning_Fine_Coarse == "Fine":
        nbins = nbins_fine
else:
        nbins = nbins_coarse
        
if mock_Type == 'cosmoSLICS':
        ID = sys.argv[1]
        if Make_Data_Vec:
                if ID !='fid':
                        print("You have set Make_Data_Vec to True, but not set the ID to fid.")
                        print("Rather than run and risk saving another cosmology as the K1000-like data vector...")
                        print("...I am just going to EXIT.")
                        sys.exit()
                Realisations = 10
                LOS_savetag = 'DatacosmoSLICS-%s-batch%s' %(Survey,Data_Vec_batch)
        else:
                Realisations = 50     # 50LOS, each PDF being the sum of the 18 regions
                LOS_savetag = 'All'
        Cov_savetag = ''      # saved as part of the covariance filename
                
if Survey == 'KiDS1000':
        #noise = ['SN0.265']
        noise = ['SN0.27', 'SN0.258', 'SN0.273', 'SN0.254', 'SN0.27']  #'SN0.28' #'NF'#
elif Survey == 'LSST':
        noise = ['SN0.28', 'SN0.28', 'SN0.28', 'SN0.28', 'SN0.28']
        
ZBcut = ['0.1-0.3', '0.3-0.5', '0.5-0.7', '0.7-0.9', '0.9-1.2']
#ZBcut = ['0.1-1.2']

DIRname_s = [] # signal directories
DIRname_n = [] # noise directories

if Read_IA:
        print(" !!! READING THE IA-CONTAMINATED MEASUREMENTS !!! ")
        IA_Tag = "IA%s_" %A_IA
else:
        IA_Tag = ''

if Calc_Combine_zbins:
        # use combinations of redshift bins as well as auto-bins
        for i in range(len(noise)):
                for j in range(i, len(noise)):
                        tmp_DIRname_s = 'MRres%s_%s100Sqdeg_%s_%s_%sGpAM_z%s_ZBcut%s' %(MRres,IA_Tag,noise[i],Mask,Survey,Survey,ZBcut[i])
                        tmp_DIRname_n = 'MRres%s_100Sqdeg_NOISE_%s_%sGpAM_z%s_ZBcut%s' %(MRres,Mask,Survey,Survey,ZBcut[i])
                        if j>i:
                                tmp_DIRname_s += '_X_ZBcut%s' %ZBcut[j]
                                # GIBLIN !!! Uncomment the following...?
                                tmp_DIRname_n += '_X_ZBcut%s' %ZBcut[j]
                                
                        if mock_Type == 'cosmoSLICS':
                                tmp_DIRname_s += '_Cosmol%s' %ID
                                
                        DIRname_s.append( tmp_DIRname_s )
                        DIRname_n.append( tmp_DIRname_n )

else:
        # just use auto bins.
        for d in range(len(noise)):
                tmp_DIRname_s = 'MRres%s_100Sqdeg_%s_%s_%sGpAM_z%s_ZBcut%s' %(MRres,noise[d],Mask,Survey,Survey,ZBcut[d])
                tmp_DIRname_n = 'MRres%s_100Sqdeg_NOISE_%s_%sGpAM_z%s_ZBcut%s' %(MRres,Mask,Survey,Survey,ZBcut[d])
                if mock_Type == 'cosmoSLICS':
                        tmp_DIRname_s += '_Cosmol%s' %ID
                DIRname_s.append( tmp_DIRname_s )
                DIRname_n.append( tmp_DIRname_n )


dir_cycle = len( DIRname_s )

if mock_Type == 'SLICS':
        
        # Read in the 217x18 matrix which indicates how the LOS&R are shuffled
        Shuffle_Config = 0     # the seed used for the shufflin'
        shffld_lr = np.load('Shuffled_SLICS-K1000-Mosaic_Configs/Shuffled_Config%s.npy' %Shuffle_Config)
        
        if Make_Data_Vec:
                Realisations = 10
                LOS_savetag = 'DataSLICS-%s-batch%s' %(Survey,Data_Vec_batch)
        else:
                # Realisations for SLICS should come out as 217
                # (los_start=74, los_end=292, no 198/199).
                Realisations = len( glob.glob( '%s/*R18*SS%s.Ekappa.npy'%(DIRname_s[0], SS[0]) ) )
                LOS_savetag = 'All-Shuffle%s' %Shuffle_Config
        ID = 'Fiducial SLICS'
        Cov_savetag = '-Shuffle%s' %Shuffle_Config      # saved as part of the covariance filename 
        
# Use the same bin_centres for every ZB/noise bin 
bin_centres = np.zeros([ len(SS), nbins ])

len_PDF_array = nbins*len(SS)*dir_cycle
Regions = 18
PDF = np.zeros([ Realisations, len_PDF_array ])
PDF_avg = np.zeros( len_PDF_array )
# curious, what do kappa (not SNR) PDFs look like?
PDF_k = np.zeros([ Regions, Realisations, len_PDF_array ])
PDF_k_avg = np.zeros([ Regions, len_PDF_array ])
use_bins_kappa = np.linspace(-0.04,0.04,nbins+1)

# also calculate and store the cumulative PDF
cumPDF = np.zeros_like( PDF )
cumPDF_avg = np.zeros_like( PDF_avg )        

index = 0      # used to append PDFs of different SS and (noise/ZBcuts) to calc overall Covariance. 
for ss in range(len(SS)):
        for d in range(dir_cycle):
                count_mask_maps = 0 # count the number of maps which have masked regions (per redshift bin)
                print( "------------- On DIRname -----------" )
                print( "------------- %s -----------" % DIRname_s[d] )

                # -------------------------------------- DEFINE SNR BINS -------------------------------------  
                if "KiDS1000" in DIRname_s[d] or 'LSST' in DIRname_s[d]:
                        # NOISY MAPS
                        if float(SS[ss])<1.5 and 'SN' in DIRname_s[d]:
                                edge = 1.5
                        elif float(SS[ss])>1.5 and float(SS[ss])<10 and 'SN' in DIRname_s[d]:
                                edge = 3.
                        elif float(SS[ss])>10 and 'SN' in DIRname_s[d]:
                                edge = 2. # (3. for original wide bins)
                                
                        # NOISELESS MAPS
                        elif float(SS[ss])>2 and float(SS[ss])<20 and 'NF' in DIRname_s[d]:
                                edge = 0.015

                        else:
                                print("KiDS1000 is specified but the binning scheme...")
                                print("...has not been defined for this smoothing scale and noise level!: SS%s" %SS[ss])
                                sys.exit()

                else:
                        print("KiDS1000 is not in the directory name: Haven't coded a binning scheme for these PDFs!")
                        sys.exit()
                
                if Binning_Fine_Coarse == "Fine":
                        # Then we're saving a very finely-binned, wide PDF
                        edge *= 4 # widen the boundaries by a factor of 4 to make sure tails
                                  # are captured by the finely binned PDF
                use_bins = np.linspace(-1.*edge, edge, nbins+1)
                        
                print("For SS %s using kappa binning: " %SS[ss])
                print(use_bins)
                bin_centres[ss,:] = use_bins[:-1] + (use_bins[1]-use_bins[0])/2.

                # -------------------------------------- FINISHED DEFINING KAPPA BINS -------------------------------------  

                # Read in Kappa maps and make a finely-binned PDF
                if Binning_Fine_Coarse=="Fine":  

                        # scroll through Mosaic regions:
                        for R in range(1,Regions+1):
                                print('------------- Region %s ---------------' %R)
                                
                                # glob the kappa maps (only used if making the finely-binned PDF)
                                filenames_s = natsorted( glob.glob( '%s/*LOS*R%s.SS%s.*Ekappa.npy'%(DIRname_s[d], R,SS[ss]) ) )
                                filenames_n = natsorted( glob.glob( '%s/*LOS*R%s.SS%s.*Ekappa.npy'%(DIRname_n[d], R,SS[ss]) ) )
                                print("The number of noise maps for Region %s is %s" %(R, len(filenames_n)) )

                                # Begin by checking if the STDEV map for this region is saved to file.
                                # If not, read in all the noise maps and make the stdev map.     
                                noise_fname = filenames_n[0].split('LOS')[0] + 'NLOS%sR%s.SS' %(len(filenames_n),R) + filenames_n[0].split('SS')[-1]
                                noise_fname = noise_fname.split('Ekappa.npy')[0] + 'NOISEstd.npy'
                                Noise_std = Read_AndOr_Make_Overall_Noise_Map(noise_fname, filenames_n)
                                print(" Region %s has a noise map of dimensionality: [%s,%s]" %(R,Noise_std.shape[0],Noise_std.shape[1]))
                                print( noise_fname )
                                
                                if Make_Data_Vec:
                                        if mock_Type == 'cosmoSLICS':
                                                # select the 0'th noise realisation for the specified batch of LOS.
                                                num_noise = int( len(filenames_s) / 50 )               # how many noise realisations have been ran
                                                filenames_s = filenames_s[::num_noise][ Data_Vec_batch*Realisations:(Data_Vec_batch+1)*Realisations ]
                                                                         # get the 0'th noise realisation per LOS
                                        else:
                                                # select first 10 LOS
                                                filenames_s = filenames_s[ Data_Vec_batch*Realisations:(Data_Vec_batch+1)*Realisations ]
                                        
                                        print( "-------------- USING THESE FILENAMES ---------------" )
                                        print( filenames_s )
                                        print( "----------------------------------------------------" )
                                        
                                if Binning_Fine_Coarse=="Fine" and len(filenames_s) != int(Realisations): # 18 for 18 regions
                                        print("We're reading in kappa maps to make the finely binned PDF...")
                                        print("...BUT there should be %s kappa maps & I only found %s! EXITING." %(Realisations, len(filenames_s)) )
                                        sys.exit()

                                # Identify masked pxls:
                                data = np.load(filenames_s[0])
                                mask_y, mask_x = np.where( data==0. )
                                Noise_std[mask_y, mask_x] = 1.          # avoid division by zero

                                if mock_Type == 'SLICS':
                                        # Shuffle the LOS & R in accordance with the shuffle matrix:
                                        filenames_s = Shuffle_filenames_LOSnR( filenames_s, R )

                                for i in range(len(filenames_s)):  # scroll through LOS
                                        #print( "Reading in file %s of %s"%(i,len(filenames_s)) )

                                        data = np.load(filenames_s[i])   # read map once & do kappa PDF
                                        #PDF_k[R-1, i, index*nbins:(index+1)*nbins ],_ = np.histogram(np.ndarray.flatten(data[data!=0.]),
                                         #                                                            use_bins_kappa,
                                         #                                                            density=False)

                                        data = np.load(filenames_s[i]) / Noise_std # read map & convert to SNR 
                                        # Set density to False, so we just get raw number of pxls, (acts as a weight). 
                                        tmp_pdf,_ = np.histogram(np.ndarray.flatten(data[data!=0.]), use_bins, density=False)
                                        # rm mask spike:
                                        tmp_pdf = Remove_Mask_Spike_From_PDF(bin_centres[ss], tmp_pdf)
                                        
                                        # save a finely-binned PDF for every realisation (LOS & R)
                                        fname_r = filenames_s[i].replace('Ekappa.npy', 'SNRPDF_%sbins.dat'%nbins)
                                        np.savetxt(fname_r, np.c_[bin_centres[ss],tmp_pdf])
                                        # sum the PDFs across regions for given LOS:
                                        PDF[i, index*nbins:(index+1)*nbins ] += tmp_pdf

                                #PDF_k_avg[R-1,:] = np.mean(PDF_k[R-1, :, index*nbins:(index+1)*nbins ], axis=0)
                        # ------------------------- FINISHED SCROLLING THROUGH R REGIONS ------------------------- #
                                        
                         # Avg the summed PDFs (where the sum is over 18 regions) across all the LOS.
                        for i in range(index*nbins, (index+1)*nbins):
                                PDF_avg[i] = np.mean( PDF[:,i] )
                                cumPDF_avg[i] = np.mean( cumPDF[:,i] )
                        # rm wee spike at centre of pdf.
                        # --- GIBLIN! CHANGED THIS; NOT NECESSARY AS SPIKE IS RM'D FROM EACH REALSN?
                        #PDF_avg[index*nbins:(index+1)*nbins] = Remove_Mask_Spike_From_PDF(bin_centres[ss],
                        #                                                                  PDF_avg[index*nbins:(index+1)*nbins] )
                        # (note that the cumulative hasn't been corrected). 
                        
                        # Assemble the filename for the saved avg PDF
                        fname_avg = fname_r.split('LOS')[0] + 'LOS%s'%(LOS_savetag) + '.SS%s.SNRPDF_%sbins.dat'%(SS[ss],nbins)
                        np.savetxt( fname_avg, np.c_[bin_centres[ss],PDF_avg[index*nbins:(index+1)*nbins] ])
                        fname_cov = fname_avg.split('LOS')[0] + 'NLOS%s%s'%(Realisations,Cov_savetag) + '.SS%s.SNRPDF_%sbins.CovMat'%(SS[ss],nbins)

                        # Also save the individual survey realisations which this avg was derived from
                        # For SLICS, this will be 217x2000 PDFs, where the LOS 've been shuffled (using seed 'Shuffle_Config')
                        # For cosmoSLICS it's 50x2000.
                        # Note this will have the same filename as the PDF_avg except for the extension (.dat-->.npy)
                        np.save( fname_avg.replace('.dat',''), PDF[:,index*nbins:(index+1)*nbins] )
                        
                else:
                # --------------------------------- COARSE PDF BINNING --------------------------- #
                        # Read in the avg finely-binned PDF and re-bin it into a smaller number of coarse bins
                        finebin_fname = glob.glob( '%s/*LOS%s.SS%s.SNRPDF_%sbins.dat'%(DIRname_s[d],LOS_savetag,
                                                                                       SS[ss],nbins_fine) )
                        if len(finebin_fname) == 0:
                                print("The finely binned PDF has not been made yet!")
                                print("Try setting Binning_Fine_Coarse to Fine, running, then back to Coarse and run again.")
                                sys.exit()
                        else:
                                print("Reading in this finely-binned PDF for re-binning:")
                                print( finebin_fname[0] )
                                
                        fine_bin_centres, fine_PDF = np.loadtxt( finebin_fname[0], usecols=(0,1), unpack=True )
                        # rebin the fine-PDF you've just read in:
                        PDF_avg[index*nbins:(index+1)*nbins], idx = Rebin_PDF(fine_bin_centres, fine_PDF, use_bins)
                        cumPDF_avg[index*nbins:(index+1)*nbins] = np.cumsum( PDF_avg[index*nbins:(index+1)*nbins] )
                        fname_avg = finebin_fname[0].replace("_%sbins"%nbins_fine, "_%sbins"%nbins)
                        

                        # Now do the same with the PDFs for the 217(50) survey realsn from SLICS(cosmoSLICS)
                        PDF_preavg = np.load(finebin_fname[0].replace('.dat','.npy'))
                        for i in range( PDF_preavg.shape[0] ): # scroll through 217(5) realsn and rebin:
                                PDF[i, index*nbins:(index+1)*nbins], idx = Rebin_PDF(fine_bin_centres, PDF_preavg[i], use_bins)
                                cumPDF[i, index*nbins:(index+1)*nbins] = np.cumsum( PDF[i, index*nbins:(index+1)*nbins] )
                        # Save the 217(50) COARSELY-BINNED survey realsn as pckld files:
                        np.save( fname_avg, PDF[:,index*nbins:(index+1)*nbins] )
                        np.save( fname_avg.replace('SNRPDF', 'SNRcumPDF'), cumPDF[:,index*nbins:(index+1)*nbins] )

                        # Following is no longer necessary - the loop above takes care of it.
                        # since we now save summed PDFs across 18 regions
                        # rather than individual LOS & R
                        '''
        		# Also read in the finely-binned PDFs for all individual realisations,
                        # for use in re-binning, and calculating the coarsely-binned covariance.
                        # Only do this if if not producing the Data vector (causes issues scrolling through all PDFs):
                        if Make_Data_Vec == False:
                                finebin_fname_all = glob.glob( '%s/*LOS*R*.SS%s.SNRPDF_%sbins.dat'%(DIRname_s[d],SS[ss],nbins_fine) )
                                finebin_fname_all = natsorted( finebin_fname_all )  # sort into order   
                                for i in range( len(finebin_fname_all) ):
                                        fine_bin_centres, fine_PDF = np.loadtxt( finebin_fname_all[i],
                                                                         usecols=(0,1), unpack=True )
                                        PDF[i, index*nbins:(index+1)*nbins], idx = Rebin_PDF(fine_bin_centres, fine_PDF, use_bins)
                                        cumPDF[i, index*nbins:(index+1)*nbins] = np.cumsum( PDF[i, index*nbins:(index+1)*nbins] )
                                        # save coarse realisation:
                                        fname_r = finebin_fname_all[i].replace("_%sbins"%nbins_fine, "_%sbins"%nbins)
                                        np.savetxt(fname_r, np.c_[bin_centres[ss],PDF[i, index*nbins:(index+1)*nbins]])
                        '''
                        # Now save the COARSELY-BINNED avg PDF (& cum PDF) you calculated above:
                        print("saving ", fname_avg)
                        np.savetxt(fname_avg,
                                   np.c_[bin_centres[ss],PDF_avg[index*nbins:(index+1)*nbins] ],
                                   header='Kappa, PDF')
                        #View_PDF( bin_centres[ss], PDF[:,index*nbins:(index+1)*nbins], PDF_avg[index*nbins:(index+1)*nbins], 0. )

                        # save the cumulative PDF
                        np.savetxt(fname_avg.replace('SNRPDF', 'SNRcumPDF'),
                                   np.c_[bin_centres[ss], cumPDF_avg[index*nbins:(index+1)*nbins] ],
                                   header='Kappa, cumPDF')

                        # Save the covariances for the individual SS and ZBcut/noise level
                        fname_cov = fname_avg.split('LOS')[0] + 'NLOS%s%s'%(Realisations,Cov_savetag) + '.SS%s.SNRPDF_%sbins.CovMat'%(SS[ss],nbins)
                        if Make_Data_Vec == False:
                                cov = np.cov(PDF[:,index*nbins:(index+1)*nbins], rowvar=False)
                                np.save(fname_cov, cov)

                                # Also save covariance for cumulative PDF
                                np.save( fname_cov.replace('SNRPDF', 'SNRcumPDF'),
                                         np.cov(cumPDF[:,index*nbins:(index+1)*nbins], rowvar=False) )

                # ------------------------- FINISHED FINE/COARSE IF STATEMENT ------------------------- #
                        
                        
                index+=1

                
# Compute the combined cov of different smoothing scales & ZBcut/noise
ss_label = 'SS'
for ss in range(len(SS)):
        ss_label += '%s-' %SS[ss]
ss_label = ss_label[:-1]
# Edit the smoothing scale tag (if only one smoothing scale is specified, this line does nothing). 
fname_cov_all = fname_cov.replace('.SS%s' %SS[-1], '.%s'%ss_label, 1)

# Now edit the redshift tag name if there's multiple bins specified.
if dir_cycle>1:
        # then we've measured the PDFs across multiple redshift bins...
        # we need to create a filename for the combined PDFs that makes sense:
        zb_label = 'zbins'
        if 'KiDS1000' in fname_cov or 'LSST' in fname_cov:
                if Calc_Combine_zbins:
                        for i in range(len(noise)):
                                # identify exactly which tomo bins are being used:
                                if ZBcut[i] == '0.1-0.3':
                                        zb1='1'
                                elif ZBcut[i] == '0.3-0.5':
                                        zb1='2'
                                elif ZBcut[i] == '0.5-0.7':
                                        zb1='3'
                                elif ZBcut[i] == '0.7-0.9':
                                        zb1='4'
                                elif ZBcut[i] == '0.9-1.2':
                                        zb1='5'
                                        
                                for j in range(i,len(noise)):
                                        if ZBcut[j] == '0.1-0.3':
                                                zb2='1'
                                        elif ZBcut[j] == '0.3-0.5':
                                                zb2='2'
                                        elif ZBcut[j] == '0.5-0.7':
                                                zb2='3'
                                        elif ZBcut[j] == '0.7-0.9':
                                                zb2='4'
                                        elif ZBcut[j] == '0.9-1.2':
                                                zb2='5'
                                                
                                        #zb_label+='%s%s-' %(i+1, j+1)
                                        zb_label+='%s%s-' %(zb1, zb2)
                                        
                else:
                
                        for dd in range(dir_cycle):
                                if ZBcut[dd] == '0.1-0.3':
                                        zb_label+='1-'
                                elif ZBcut[dd] == '0.3-0.5':
                                        zb_label+='2-'
                                elif ZBcut[dd] == '0.5-0.7':
                                        zb_label+='3-'
                                elif ZBcut[dd] == '0.7-0.9':
                                        zb_label+='4-'
                                elif ZBcut[dd] == '0.9-1.2':
                                        zb_label+='5-'
                zb_label = zb_label[:-1]
                zb_label += '.'
        else:
                print("I haven't coded up what the multiple-zbin covariance savename should be for anything other than KiDS1000-like mocks. This is not the KiDS1000 mocks, but you have set multiple redshifts. Edit this bit of code to continue.")
                sys.exit()
                
        # Inseting zb_label into the filesavename if there's multiple redshifts
        # Saving the cov for multiple zbins in the directory specified by the final ZBcut flag.
        fname_cov_all = fname_cov_all.split('bins.')[0] + 'bins.%sCovMat' %zb_label + fname_cov_all.split('CovMat')[-1]

cov = np.cov(PDF, rowvar=False)

if Binning_Fine_Coarse == "Coarse" and Make_Data_Vec == False:
        np.save(fname_cov_all, cov)
        np.save(fname_cov_all.replace('SNRPDF', 'SNRcumPDF'), np.cov(cumPDF) )

        # Save all PDFs as a pickled data file, to be read in to calc.
        # the combined cov of xi+/- and PDFs
        np.save( fname_cov_all.replace('.CovMat', ''), PDF )
	
