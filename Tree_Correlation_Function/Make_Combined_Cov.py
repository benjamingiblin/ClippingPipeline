# 06/11/2019, B. M. Giblin, Postdoc, Edinburgh
# Read in the clipped & unclipped xi+ for all redshift bins and lines of sight, calculate combined covariance.

import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import sys
from astropy.io import fits

def View_xi(theta, xi_c_all, xi_c_avg, xi_uc_all, xi_uc_avg, pm):
        plt.figure()
        for i in range(xi_c_all.shape[0]):
                plt.plot(theta, theta*1e4*xi_c_all[i,:], color='dimgrey')
        for i in range(xi_uc_all.shape[0]):
                plt.plot(theta, theta*1e4*xi_uc_all[i,:], color='yellow')
        plt.plot(theta, theta*1e4*xi_c_avg, color='black', linewidth=3, label='Clipped')
        plt.plot(theta, theta*1e4*xi_uc_avg, color='orange', linewidth=3, label='Unclipped')

        plt.xscale('log')
        plt.xlabel(r'$\theta$ [arcmin]')
        plt.ylabel(r'$\theta \times \xi_%s \, [10^{-4} \rm{arcmin}]$' %pm)
        #plt.legend(loc='best', frameon=False)
        plt.show()
        return

mock_Type = 'SLICS'    # 'cosmoSLICS' 'SLICS'
Auto_Or_Cross = 'Cross' #'Cross' # If set to 'Auto', only reads from the 5 auto-correlation bins
                       # If set to 'Cross', reads from auto AND cross bins (15 combinations)

if mock_Type == 'cosmoSLICS':
        ID = sys.argv[1]
        cosmol_flag='_Cosmol%s' %ID
        Realisations_c=50
        Realisations_uc=50
else:
        cosmol_flag = ''

# KiDS1000-like specs
Survey = 'KiDS1000'
MRres = 'MRres140.64arcs_' #'60arcs_'
Mask = 'Mosaic'
#noise = [ 'SN0.27', 'SN0.258', 'SN0.273', 'SN0.254', 'SN0.27']  #['SN0.265'] 
#ZBcut = ['0.1-0.3', '0.3-0.5', '0.5-0.7', '0.7-0.9', '0.9-1.2'] #['0.1-1.2']
# use this for ZBcut cuts:
noise = [ 'SN0.254', 'SN0.27']
ZBcut = ['0.7-0.9', '0.9-1.2']

NLOS_flag = 3906 #715

# My LSST-like specs
#Survey = 'LSST'
#MRres = 'MRres60arcs_'
#noise = [ 'SN0.28', 'SN0.28', 'SN0.28', 'SN0.28', 'SN0.28']
#ZBcut = ['0.1-0.3', '0.3-0.5', '0.5-0.7', '0.7-0.9', '0.9-1.2']
#NLOS_flag = 616

# Chris Davies' LSST-like specs
#Survey = 'LSST'
#MRres = 'MRres10arcs_'
#noise = ['SN0.28']  
#ZBcut = ['0.6-1.4']

if Mask == "Mosaic":
        ftag = 'Mosaic'
elif Mask == 'NoMask':
        ftag = 'test'

# Variable for calculating the combined xi+/- AND PDF
Combine_with_PDF = True
PDF_SS_label = 'SS1.408-2.816-5.631' 
PDF_nbins = 4
PDF_DIR = '/home/bengib/Clipping_Pipeline/Mass_Recon/%s100Sqdeg_%s_%s_%sGpAM_z%s_ZBcut%s/' %(MRres,noise[-1],Mask,
                                                                                            Survey, Survey,
                                                                                            ZBcut[-1])
LOS_savetag = 'All'

# Assemble the directories to cycle through
HOME_DIR = '/home/bengib/Clipping_Pipeline/Tree_Correlation_Function'
DIRname = []

if Auto_Or_Cross == 'Auto':
        for d in range(len(ZBcut)):
                DIRname.append('%s/%s100Sqdeg_%s_%s_%sGpAM_z%s_ZBcut%s%s/ThBins9' %(HOME_DIR, MRres, noise[d],
                                                                                    Mask, Survey, Survey,
                                                                                    ZBcut[d], cosmol_flag))

elif Auto_Or_Cross == 'Cross':
        for i in range(len(ZBcut)):
                for j in range(i, len(ZBcut)):
                        DIRname.append('%s/%s100Sqdeg_%s_%s_%sGpAM_z%s_ZBcut%s_X_ZBcut%s%s/ThBins9' %(HOME_DIR, MRres,
                                                                                                      noise[i],Mask,
                                                                                                      Survey,
                                                                                                      Survey,
                                                                                                      ZBcut[i],
                                                                                                      ZBcut[j],
                                                                                                      cosmol_flag))
                
dir_cycle = len(DIRname)
SS = [2.816] #1, 9.33, 84.85 
        
if mock_Type == 'SLICS':
        filenames_c = sorted( glob.glob( '%s/%s_%s.%sGpAM.LOS*.SS%s.rCLIP_X3sigma.CorrFun.asc'%(DIRname[0],noise[0],
                                                                                                ftag,Survey,
                                                                                                SS[0]) ) )
        Realisations_c = len(filenames_c)
        filenames_uc = sorted( glob.glob( '%s/%s_%s.%sGpAM.LOS*.ORIG.CorrFun.asc'%(DIRname[0],noise[0],
                                                                                   ftag, Survey)) )
        Realisations_uc = len(filenames_uc)
        
        ID = 'Fiducial SLICS'
        
nbins=9
len_array = nbins*len(SS)*dir_cycle
# clipped
xi_p_c_all = np.zeros([ Realisations_c, len_array ])   # xi+
xi_m_c_all = np.zeros([ Realisations_c, len_array ])   # xi-
xi_p_c_avg = np.zeros( len_array )
xi_m_c_avg = np.zeros( len_array )
# unclipped
xi_p_uc_all = np.zeros([ Realisations_uc, len_array ]) # xi+
xi_m_uc_all = np.zeros([ Realisations_uc, len_array ]) # xi-
xi_p_uc_avg = np.zeros_like( xi_p_c_avg )
xi_m_uc_avg = np.zeros_like( xi_p_c_avg )


index = 0      # used to append PDFs of different SS and (noise/ZBcuts) to calc overall Covariance. 
for ss in range(len(SS)):
        for d in range(dir_cycle):
                print( "------------- On mock-type %s, Cosmol %s, SS %s, DIRname %s -----------" %(mock_Type, ID,SS[ss],DIRname[d]) )
                                                                                            
                filenames_c = sorted( glob.glob( '%s/*_%s.%sGpAM.LOS*.SS%s.rCLIP_X3sigma.CorrFun.asc'%(DIRname[d],
                                                                                                       ftag,Survey,
                                                                                                       SS[ss])) )
                filenames_uc = sorted( glob.glob( '%s/*_%s.%sGpAM.LOS*.ORIG.CorrFun.asc'%(DIRname[d],ftag,Survey)) )

                # if statements to exit if number of files across zbins differs
                if len(filenames_c) != Realisations_c:
                        print("There's %s clipped filenames in first DIRname." %Realisations_c)
                        print("But only found %s in DIRname number %s! EXITING." %(len(filenames_c), d) )
                        sys.exit()
                if len(filenames_uc) != Realisations_uc:
                        print("There's %s unclipped filenames in first DIRname." %Realisations_uc)
                        print("But only found %s in DIRname number %s! EXITING." %(len(filenames_uc), d) )
                        sys.exit()
                        
                
                for i in range(len(filenames_c)):
                        print( "Reading in clipped file %s of %s"%(i,len(filenames_c)) )
                        theta, xi_p_c_all[i, index*nbins:(index+1)*nbins] = np.loadtxt( filenames_c[i], usecols=(1,3), unpack=True)
                        theta, xi_m_c_all[i, index*nbins:(index+1)*nbins] = np.loadtxt( filenames_c[i], usecols=(1,4), unpack=True)
                        
                for i in range(len(filenames_uc)):
                        print( "Reading in unclipped file %s of %s"%(i,len(filenames_uc)) )
                        xi_p_uc_all[i, index*nbins:(index+1)*nbins] = np.loadtxt( filenames_uc[i], usecols=(3,), unpack=True)
                        xi_m_uc_all[i, index*nbins:(index+1)*nbins] = np.loadtxt( filenames_uc[i], usecols=(4,), unpack=True)

                for i in range(index*nbins, (index+1)*nbins):
                        xi_p_c_avg[i] =  np.mean( xi_p_c_all[:,i] )
                        xi_p_uc_avg[i] = np.mean( xi_p_uc_all[:,i] )
                        xi_m_c_avg[i] =  np.mean( xi_m_c_all[:,i] )
                        xi_m_uc_avg[i] = np.mean( xi_m_uc_all[:,i] )
                        
                index+=1

                
# Plot the xi+/-
#for i in range(dir_cycle):
#        View_xi( theta,
#                  xi_p_c_all[:,i*nbins:(i+1)*nbins], xi_p_c_avg[i*nbins:(i+1)*nbins],
#                  xi_p_uc_all[:,i*nbins:(i+1)*nbins], xi_p_uc_avg[i*nbins:(i+1)*nbins], '+' )
#        View_xi( theta,
#                  xi_m_c_all[:,i*nbins:(i+1)*nbins], xi_m_c_avg[i*nbins:(i+1)*nbins],
#                  xi_m_uc_all[:,i*nbins:(i+1)*nbins], xi_m_uc_avg[i*nbins:(i+1)*nbins], '-' )
        

# Prepare the filename for the combined cov of different smoothing scales, ZBcut/noise, clipped & unclipped:
OUTDIR_c = '%s/NLOS%s' %(DIRname[0],len(filenames_c))
OUTDIR_uc = '%s/NLOS%s' %(DIRname[0],len(filenames_uc))
# make output directories if they dont exist
if not os.path.exists(OUTDIR_c):
    os.makedirs(OUTDIR_c)
if not os.path.exists(OUTDIR_uc):
    os.makedirs(OUTDIR_uc)
                         
fname_cov_p_c = '%s/%s_%s.%sGpAM.NLOS%s.SS%s.rCLIP_X3sigma.CxC+.CovMat.npy' %(OUTDIR_c, noise[0], ftag, Survey, len(filenames_c),SS[0])
fname_cov_p_uc = '%s/%s_%s.%sGpAM.NLOS%s.ORIG.UCxUC+.CovMat.npy' %(OUTDIR_uc, noise[0], ftag, Survey, len(filenames_uc))
fname_cov_m_c = '%s/%s_%s.%sGpAM.NLOS%s.SS%s.rCLIP_X3sigma.CxC-.CovMat.npy' %(OUTDIR_c, noise[0], ftag, Survey, len(filenames_c),SS[0])
fname_cov_m_uc = '%s/%s_%s.%sGpAM.NLOS%s.ORIG.UCxUC-.CovMat.npy' %(OUTDIR_uc, noise[0], ftag, Survey, len(filenames_uc))

ss_label = 'SS'
for ss in range(len(SS)):
        ss_label += '%s-' %SS[ss]
ss_label = ss_label[:-1]

# Edit the smoothing scale tag (if only one smoothing scale is specified, this line does nothing). 
fname_cov_p_c = fname_cov_p_c.replace('.SS%s' %SS[-1], '.%s'%ss_label, 1)
fname_cov_m_c = fname_cov_m_c.replace('.SS%s' %SS[-1], '.%s'%ss_label, 1)

# Now edit the redshift tag name if there's multiple bins specified.
if dir_cycle>1:
        # then we've measured the PDFs across multiple redshift bins...
        # we need to create a filename for the combined PDFs that makes sense:
        zb_label = 'zbins'
        if 'KiDS1000' in DIRname[0] or 'LSST' in DIRname[0]:
                if Auto_Or_Cross == 'Auto':
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

                elif Auto_Or_Cross == 'Cross':
                        for i in range(len(ZBcut)):
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
                                for j in range(i, len(ZBcut)):
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

                                        zb_label+='%s%s_' %(zb1, zb2)
                                        #zb_label += '%s%s_' %(i+1,j+1)
                                        
                zb_label = zb_label[:-1]
                # Injecting zb_label into the filesavename if there's multiple redshifts
                # Saving the cov for multiple zbins in the directory specified by the final ZBcut flag.
                fname_cov_p_c = fname_cov_p_c.split('.npy')[0] + '.%s.npy' %zb_label
                fname_cov_p_uc = fname_cov_p_uc.split('.npy')[0] + '.%s.npy' %zb_label
                fname_cov_m_c = fname_cov_m_c.split('.npy')[0] + '.%s.npy' %zb_label
                fname_cov_m_uc = fname_cov_m_uc.split('.npy')[0] + '.%s.npy' %zb_label
                
        else:
                print("I haven't coded up what the multiple-zbin covariance savename should be for anything other than KiDS1000-like and LSST-like mocks. This is neither, but you have set multiple redshifts. Edit this bit of code to continue.")
                sys.exit()

cov_p_c = np.cov(xi_p_c_all, rowvar=False)
cov_p_uc = np.cov(xi_p_uc_all, rowvar=False)
cov_m_c = np.cov(xi_m_c_all, rowvar=False)
cov_m_uc = np.cov(xi_m_uc_all, rowvar=False)

np.save(fname_cov_p_c, cov_p_c)
np.save(fname_cov_p_uc, cov_p_uc)
np.save(fname_cov_m_c, cov_m_c)
np.save(fname_cov_m_uc, cov_m_uc)

# Combine the xi+ and xi-
xi_c_all = np.concatenate((xi_p_c_all,xi_m_c_all), axis=1)
xi_uc_all = np.concatenate((xi_p_uc_all,xi_m_uc_all), axis=1)
# calc xi+/- cov
cov_c = np.cov( xi_c_all, rowvar=False )
cov_uc = np.cov( xi_uc_all, rowvar=False )
# prepare savename
fname_cov_c = fname_cov_p_c.replace('CxC+', 'CxC+-', 1)
fname_cov_uc = fname_cov_p_uc.replace('UCxUC+', 'UCxUC+-', 1)
# save cov's
np.save(fname_cov_c, cov_c)
np.save(fname_cov_uc, cov_uc)

# Compute combined covariance for clipped and unclipped:
if Realisations_c == Realisations_uc:
        # equal number of clipped and unclipped files. Safe to calc combined covariance.

        # xi+ clipped and unclipped
        xi_p_cuc_all = np.concatenate((xi_p_c_all,xi_p_uc_all), axis=1)
        cov_p_cuc = np.cov( xi_p_cuc_all, rowvar=False ) 
        fname_cov_p_cuc = fname_cov_p_c.replace('CxC+', 'ClipAndUnclip+', 1)
        np.save(fname_cov_p_cuc, cov_p_cuc)

        # xi- clipped and unclipped
        xi_m_cuc_all = np.concatenate((xi_m_c_all,xi_m_uc_all), axis=1)
        cov_m_cuc = np.cov( xi_m_cuc_all, rowvar=False )
        fname_cov_m_cuc = fname_cov_m_c.replace('CxC-', 'ClipAndUnclip-', 1)
        np.save(fname_cov_m_cuc, cov_m_cuc)

        # xi+ AND xi- clipped and unclipped
        xi_cuc_all = np.concatenate((xi_c_all,xi_uc_all), axis=1)
        cov_cuc = np.cov( xi_cuc_all, rowvar=False )
        fname_cov_cuc = fname_cov_m_c.replace('CxC-', 'ClipAndUnclip+-', 1)
        np.save(fname_cov_cuc, cov_cuc)
        
else:
        print("Number of clipped and unclipped files is different. Wont calc. combined cov of them both. EXITING.")
        sys.exit()



PDF_fname = PDF_DIR + '%s_%s.%sGpAM.NLOS%s.%s.SNRPDF_%sbins.%s.npy' %(noise[-1],
                                                                      ftag,Survey,
                                                                      NLOS_flag,
                                                                      PDF_SS_label,
                                                                      PDF_nbins,zb_label) 
if Combine_with_PDF:
        PDFs = np.load( PDF_fname )
        # Combined PDF with UNCLIPPED xi_+/-
        PDF_xi_uc_all = np.concatenate((PDFs,xi_uc_all), axis=1)
        cov_PDF_uc = np.cov( PDF_xi_uc_all, rowvar=False )
        fname_cov_PDF_uc = fname_cov_uc.replace('UCxUC+-', 'PDF-6.6arcmin-UCxUC+-', 1)
        np.save( fname_cov_PDF_uc, cov_PDF_uc )



        
def Convert_Cov_2_CCC(cov):
        # cross-corr coeff matrix
        ccc =  np.empty_like(cov)
        for i in range(cov.shape[0]):
                for j in range(cov.shape[1]):
                        ccc[i,j] = cov[i,j] / np.sqrt(  abs(cov[i,i]*cov[j,j]) )
        return ccc
                        
def Plot_CCC(ccc, theta, title, savename):
        plt.figure()
        plt.imshow(ccc, vmin=-1, vmax=1., origin='lower', interpolation='nearest') #,
                   #extent = [theta.min(), theta.max(), theta.min(), theta.max()])
        #plt.ylabel(r'$\theta$ [arcmin]')
        #plt.xlabel(r'$\theta^{\prime}$ [arcmin]')
        plt.colorbar()
        plt.title(title)
        #plt.savefig('%sCCC_Mat.png'%(savename))
        plt.show()                                                                                                        
        return

# xi+
#Plot_CCC( Convert_Cov_2_CCC(cov_p_c), theta, r'$\xi_+^{\rm clip}$', None )
#Plot_CCC( Convert_Cov_2_CCC(cov_p_uc), theta, r'$\xi_+^{\rm unclip}$', None )
#Plot_CCC( Convert_Cov_2_CCC(cov_p_cuc), theta, r'$\xi_+^{\rm clip \, and \, unclip}$', None )

# xi+ AND xi-
#Plot_CCC( Convert_Cov_2_CCC(cov_c), theta, r'$\xi_\pm^{\rm clip}$', None )
#Plot_CCC( Convert_Cov_2_CCC(cov_uc), theta, r'$\xi_\pm^{\rm unclip}$', None )
#Plot_CCC( Convert_Cov_2_CCC(cov_cuc), theta, r'$\xi_\pm^{\rm clip \, and \, unclip}$', None )

# PDF AND xi+^unclip AND xi-^unclip
#Plot_CCC( Convert_Cov_2_CCC(cov_PDF_uc), theta, r'PDF & $\xi_\pm^{\rm unclip}$', None )  
