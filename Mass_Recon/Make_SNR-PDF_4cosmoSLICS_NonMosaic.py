# 06/11/2019, B. M. Giblin, Postdoc, Edinburgh
# Read in the Kappa Maps for the cosmoSLICS cosmologies, convert to PDFs and save.
# THIS ONE WORKS ON THE NON-MOSAIC MOCKS (doesn't scroll through 18regions)

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
        new_pdf = tmp_interp( kappa ) # same as old for all idx except central 4 idx.
        
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


# ------------------------------------------------------------------------ #

# Define binning for the PDF
# If "Fine" it reads in a map and makes a finely-binned PDF
# If "Coarse" it reads in the finely-binned PDF and re-bins it,
# to make a coarsely-binned PDF.
# This saves us storing all the kappa maps (several TB).                                                                                
Binning_Fine_Coarse = "Coarse" #"Coarse"

# ------------------------------------------------------------------------ #   

mock_Type = 'cosmoSLICS'    # 'cosmoSLICS' 'SLICS'
Survey = 'KiDS1000'    # 'LSST' 'KiDS1000'
Mask = 'NoMask'        # NoMask'
Read_IA = False        # if True read the measurements for the (fid-cosmol) IA-contam. measurements
Read_SLC = True       # read SLC (IA=0 with SN corr'd with delta); CANNOT HAVE BOTH _IA/_SLC True AT ONCE!
if Read_IA or Read_SLC:
        A_IA = sys.argv[2] #"0.0"           # IA/SLC amplitude

Calc_Combine_zbins = True   # If True, will read in the PDFs for combinations of redshift bins

Make_Data_Vec = False   # If True, 
                       # it will just read & avg the 1st 10 maps & avg them
                       # to make something with the noise levels of KiDS1000

Scale_To_Mosaic = True # Scale the pxl counts so as to be the same area as the mosaic mocks
if Scale_To_Mosaic:
        npxl_mosaic = float( np.loadtxt('/home/bengib/KiDS1000_Data/Masks_SLICS_Regions/MasksArea_PxlCount_TotArea_140.64arcs.txt', 
                                        usecols=(0,), unpack=True) )


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
                Realisations = 50     # 50LOS 
                LOS_savetag = 'All'

                
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
elif Read_SLC:
        print(" !!! READING THE SLC-CONTAMINATED MEASUREMENTS !!! ")
        IA_Tag = "SLC%s_" %A_IA # It's called IA_Tag, but its really SLC
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
        if Make_Data_Vec:
                Realisations = 10
                LOS_savetag = 'DataSLICS-%s-batch%s' %(Survey,Data_Vec_batch)
        else:
                Realisations = len( glob.glob( '%s/*SS%s.Ekappa.npy'%(DIRname_s[0], SS[0]) ) )
                LOS_savetag = 'All'
        ID = 'Fiducial SLICS'

# Use the same bin_centres for every ZB/noise bin 
bin_centres = np.zeros([ len(SS), nbins ])

len_PDF_array = nbins*len(SS)*dir_cycle
PDF = np.zeros([ Realisations, len_PDF_array ])
PDF_avg = np.zeros( len_PDF_array )

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
                        
                        # glob the kappa maps (only used if making the finely-binned PDF)
                        filenames_s = natsorted( glob.glob( '%s/*LOS*.SS%s.*Ekappa.npy'%(DIRname_s[d], SS[ss]) ) )
                        filenames_n = natsorted( glob.glob( '%s/*LOS*.SS%s.*Ekappa.npy'%(DIRname_n[d], SS[ss]) ) )
                        print("The number of noise maps is %s" %len(filenames_n) )

                        # Begin by checking if the STDEV map for this region is saved to file.
                        # If not, read in all the noise maps and make the stdev map.     
                        noise_fname = filenames_n[0].split('LOS')[0] + 'NLOS%s.SS' %(len(filenames_n)) + filenames_n[0].split('SS')[-1]
                        noise_fname = noise_fname.split('Ekappa.npy')[0] + 'NOISEstd.npy'
                        Noise_std = Read_AndOr_Make_Overall_Noise_Map(noise_fname, filenames_n)
                        
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
                                
                        if Binning_Fine_Coarse=="Fine" and len(filenames_s) != Realisations: # 18 for 18 regions
                                print("We're reading in kappa maps to make the finely binned PDF...")
                                print("...BUT there should be %s kappa maps & I only found %s! EXITING." %(Realisations, len(filenames_s)) )
                                sys.exit()


                        # Identify masked pxls & get num. pxls in map
                        data = np.load(filenames_s[0])
                        npxl_map = data.shape[0] * data.shape[1]
                        #print("There are %s empty pxls in kappa map of %s pxls" %(len(np.where( data==0. )[0]), npxl_map) )
                        #mask_y, mask_x = np.where( data==0. )
                        #Noise_std[mask_y, mask_x] = 1.          # avoid division by zero                                

                        for i in range(len(filenames_s)):
                                #print( "Reading in file %s of %s"%(i,len(filenames_s)) )
                                data = np.load(filenames_s[i]) / Noise_std  # read map & convert to SNR

                                # Set density to False, so the num_pxls (weight) is accounted for when we avg.
                                tmp_pdf,_ = np.histogram(np.ndarray.flatten(data[data!=0.]), use_bins, density=False)
                                if Scale_To_Mosaic:
                                        tmp_pdf = tmp_pdf*(npxl_mosaic/npxl_map) # scale num pxls in map up to num. in K1000 mosaic. 
                                # store:
                                PDF[i, index*nbins:(index+1)*nbins] = tmp_pdf
                                cumPDF[i, index*nbins:(index+1)*nbins] = np.cumsum( tmp_pdf )

                                # save a finely-binned PDF for every realisation
                                fname_r = filenames_s[i].replace('Ekappa.npy', 'SNRPDF_%sbins.dat'%nbins)
                                np.savetxt(fname_r, np.c_[bin_centres[ss],PDF[i, index*nbins:(index+1)*nbins]])

                                        
                         # Avg the summed PDFs (weighted avg given num_pxls is accounted for in the sum)
                        for i in range(index*nbins, (index+1)*nbins):
                                PDF_avg[i] = np.mean( PDF[:,i] )
                                cumPDF_avg[i] = np.mean( cumPDF[:,i] )
                        # rm wee spike at centre of pdf.
                        PDF_avg[index*nbins:(index+1)*nbins] = Remove_Mask_Spike_From_PDF(bin_centres[ss],
                                                                                          PDF_avg[index*nbins:(index+1)*nbins] )
                        # (note that the cumulative hasn't been corrected). 
                        
                        # Assemble the filename for the saved avg PDF
                        fname_avg = fname_r.split('LOS')[0] + 'LOS%s'%(LOS_savetag) + '.SS%s.SNRPDF_%sbins.dat'%(SS[ss],nbins)
                        np.savetxt( fname_avg, np.c_[bin_centres[ss],PDF_avg[index*nbins:(index+1)*nbins] ])
                        fname_cov = fname_avg.split('LOS')[0] + 'NLOS%s'%Realisations + '.SS%s.SNRPDF_%sbins.CovMat'%(SS[ss],nbins)

                        if Read_IA:
                                # Create an IA-contaminated mosaic data vector:
                                # A_IA=0 data:                                                                            
                                fname_ia0 = fname_avg.replace(IA_Tag, 'IA0.0_')
                                _, pdf_ia0 = np.loadtxt(fname_ia0, usecols=(0,1), unpack=True)
                                # bias the data vector      
                                ia_bias = PDF_avg[index*nbins:(index+1)*nbins] - pdf_ia0 # additive
                                ia_bias2 = PDF_avg[index*nbins:(index+1)*nbins]/pdf_ia0  # scalar bias
                                ia_bias2[np.isfinite(ia_bias2)==False]=1.                # kill places with division by 0.      

                                # save the pure bias vector (for sys-emul in MCMC to train on)
                                # just use the additive, multplicative is noisy.                           
                                fname_purebias = fname_avg.replace('.dat', '.biasIA%s-additive.dat' %A_IA)
                                np.savetxt( fname_purebias, np.c_[bin_centres[ss], ia_bias])
                                print("saved this pure IA data vector:")
                                print(fname_purebias)
                                print("---------------")

                                # take the diff of IA0 and non-mosaic CS to get the bias by
                                # source-lens clustering (SLC):
                                if A_IA == '0.0':
                                        fname_noSLC = fname_ia0.replace('IA0.0_', '')
                                        _, pdf_noSLC = np.loadtxt(fname_noSLC, usecols=(0,1), unpack=True)
                                        slc_bias = pdf_ia0 - pdf_noSLC # additive
                                        slc_bias2= pdf_ia0 / pdf_noSLC # scalar (but not currently saving this)
                                        slc_bias2[np.where(np.isfinite(slc_bias2) == False)] = 1. # deal with bad pxls:
                                        slc_bias2[np.where(np.isfinite(slc_bias2)>1e7)] = 1.
                                
                                        # save the SLC biases
                                        fname_slcbias = fname_purebias.replace('.biasIA0.0-add','.biasSLC-add')
                                        np.savetxt( fname_slcbias, np.c_[bin_centres[ss], slc_bias])
                                        print("saved this pure SLC data vector:")
                                        print(fname_purebias)
                                        print("---------------")
                                        
                                # Do 2 contaminations: one SLICS, one cosmoSLICS:
                                for cs_or_s in ("cosmoSLICS", "SLICS"):
                                        print("Contaminating %s with %s" %(cs_or_s,IA_Tag))
                                        # data to contaminate:
                                        fname_msc = fname_avg.replace(IA_Tag, '').replace(Mask, 'Mosaic').replace('_test', '_Mosaic')
                                        if cs_or_s=="SLICS":
                                                fname_msc = fname_msc.replace('_Cosmolfid','').replace('LOSAll','LOSAll-Shuffle0')
                                        bins_msc, pdf_msc = np.loadtxt(fname_msc, usecols=(0,1), unpack=True)

                                        # save the mosaic pdf bias'd by the IAs:
                                        fname_bias = fname_msc.replace('.dat', '.biasIA%s-added.dat' %A_IA)
                                        fname_bias2 = fname_msc.replace('.dat', '.biasIA%s-scaled.dat' %A_IA)
                                        np.savetxt( fname_bias, np.c_[bin_centres[ss], pdf_msc+ia_bias])   # additive bias
                                        np.savetxt( fname_bias2, np.c_[bin_centres[ss], pdf_msc*ia_bias2]) # scaled bias
                                        '''
                                        print(" Saved these IA-contaminated mosaic data vector: ")
                                        print(fname_bias)
                                        print(fname_bias2)
                                        print("---------------")
                                        '''
                                        if A_IA	== '0.0':
                                                # and now contaminated by SLC:
                                                fname_slc = fname_msc.replace('.dat', '.biasSLC-added.dat')
                                                fname_slc2 = fname_msc.replace('.dat', '.biasSLC-scaled.dat')
                                                # additive/scalar:
                                                np.savetxt( fname_slc, np.c_[bin_centres[ss], pdf_msc+slc_bias])
                                                np.savetxt( fname_slc2, np.c_[bin_centres[ss], pdf_msc*slc_bias2]) 
                                                print(" Saved these SLC-contaminated mosaic data vector: ")
                                                print(fname_slc)
                                                print(fname_slc2)
                                        
                else:
                # --------------------------------- COARSE PDF BINNING --------------------------- #
                        # Read in the avg finely-binned PDF and re-bin it into a smaller number of coarse bins
                        finebin_fname = glob.glob( '%s/*LOSAll.SS%s.SNRPDF_%sbins.dat'%(DIRname_s[d],SS[ss],nbins_fine) )
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

        		# Also read in the finely-binned PDFs for all individual realisations,
                        # for use in re-binning, and calculating the coarsely-binned covariance.
                        # Only do this if if not producing the Data vector (causes issues scrolling through all PDFs):
                        if Make_Data_Vec == False:
                                finebin_fname_all = glob.glob( '%s/*LOS*.SS%s.SNRPDF_%sbins.dat'%(DIRname_s[d],SS[ss],nbins_fine) )
                                # think these lines are necessary?
                                finebin_fname_all.remove(finebin_fname[0])             # Get rid of the avg
                                finebin_fname_all = [ x for x in finebin_fname_all if "batch" not in x ] # get rid of Data measurement
                                                                                                         # (works for fid cosmol)
                                finebin_fname_all = natsorted( finebin_fname_all )  # sort into order   
                                for i in range( len(finebin_fname_all) ):
                                        fine_bin_centres, fine_PDF = np.loadtxt( finebin_fname_all[i],
                                                                         usecols=(0,1), unpack=True )
                                        PDF[i, index*nbins:(index+1)*nbins], idx = Rebin_PDF(fine_bin_centres, fine_PDF, use_bins)
                                        cumPDF[i, index*nbins:(index+1)*nbins] = np.cumsum( PDF[i, index*nbins:(index+1)*nbins] )
                                        # save coarse realisation:
                                        fname_r = finebin_fname_all[i].replace("_%sbins"%nbins_fine, "_%sbins"%nbins)
                                        np.savetxt(fname_r, np.c_[bin_centres[ss],PDF[i, index*nbins:(index+1)*nbins]])

                        print("saving ", fname_avg)
                        np.savetxt(fname_avg,
                                   np.c_[bin_centres[ss],PDF_avg[index*nbins:(index+1)*nbins] ],
                                   header='Kappa, PDF')
                        #View_PDF( bin_centres[ss], PDF[:,index*nbins:(index+1)*nbins], PDF_avg[index*nbins:(index+1)*nbins], 0. )

                        # save the cumulative PDF
                        np.savetxt(fname_avg.replace('SNRPDF', 'SNRcumPDF'),
                                   np.c_[bin_centres[ss], cumPDF_avg[index*nbins:(index+1)*nbins] ],
                                   header='Kappa, cumPDF')

                        if Read_IA:
                                # Rebin the pure bias vector (for sys-emul in MCMC to train on) for both IA & SLC(IA=0 only)
                                #if A_IA == '0.0':
                                #        loop=['IA%s'%A_IA,'SLC'] # there's only SLC bias for A_IA=0 directories
                                #else:
                                #        loop=['IA%s'%A_IA]
                                loop=['IA%s'%A_IA] # no longer saving a pure SLC bias; this gets done in Plot_*.py
                                        
                                for ibias in loop: 
                                        # get the filename:
                                        fname_purebias = fname_avg.replace('.dat', '.bias%s-additive.dat' %ibias) # get fname
                                        fname_purebias = fname_purebias.replace("_%sbins"%nbins, "_%sbins"%nbins_fine)
                                        # read it in & rebin:
                                        _, purebias = np.loadtxt(fname_purebias, usecols=(0,1), unpack=True) 
                                        purebias_coarse,_ = Rebin_PDF(fine_bin_centres, purebias, use_bins)
                                        # make new fname & save:
                                        fname_purebias = fname_purebias.replace("_%sbins"%nbins_fine, "_%sbins"%nbins) 
                                        np.savetxt( fname_purebias, np.c_[bin_centres[ss], purebias_coarse]) 
                                        print("saved this pure IA data vector:")
                                        print(fname_purebias)
                                        print("---------------")
                                
                                        # Re-bin the IA-contaminated mosaic data vector:
                                        # Do this twice: one for cosmoSLICS and one for SLICS:
                                        for cs_or_s in ("cosmoSLICS", "SLICS"):
                                                print("Contaminating %s with %s" %(cs_or_s,ibias))
                                                fname_bias_tmp = fname_avg.replace(IA_Tag, '').replace(Mask,
                                                                                                       'Mosaic').replace('_test', '_Mosaic')
                                                if cs_or_s=="SLICS":
                                                        fname_bias_tmp = fname_bias_tmp.replace('_Cosmolfid',
                                                                                                '').replace('LOSAll','LOSAll-Shuffle0')
                                                fname_bias = fname_bias_tmp.replace("_%sbins"%nbins, "_%sbins"%nbins_fine).replace('.dat', '.bias%s-added.dat' %ibias)
                                                fname_bias2 = fname_bias_tmp.replace("_%sbins"%nbins, "_%sbins"%nbins_fine).replace('.dat', '.bias%s-scaled.dat' %ibias)
                                
                                                _, pdf_bias =np.loadtxt(fname_bias, usecols=(0,1), unpack=True)  # additive bias
                                                _, pdf_bias2 =np.loadtxt(fname_bias2, usecols=(0,1), unpack=True) # scaled bias

                                                # Re-bin & save the IA-contaminated:
                                                pdf_bias_coarse,_ = Rebin_PDF(fine_bin_centres, pdf_bias, use_bins)
                                                pdf_bias_coarse2,_ = Rebin_PDF(fine_bin_centres, pdf_bias2, use_bins)
                                
                                                fname_bias = fname_bias.replace("_%sbins"%nbins_fine, "_%sbins"%nbins)
                                                fname_bias2 = fname_bias2.replace("_%sbins"%nbins_fine, "_%sbins"%nbins)
                                
                                                np.savetxt( fname_bias, np.c_[bin_centres[ss], pdf_bias_coarse])
                                                np.savetxt( fname_bias2, np.c_[bin_centres[ss], pdf_bias_coarse2])
                                                print(" Saved this contaminated mosaic data vector: ")
                                                print(fname_bias)
                                                print(fname_bias2)
                                        print("---------------")
                        

                        # Save the covariances for the individual SS and ZBcut/noise level
                        fname_cov = fname_avg.split('LOS')[0] + 'NLOS%s'%Realisations + '.SS%s.SNRPDF_%sbins.CovMat'%(SS[ss],nbins)
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
	
