# 03/08/2021, B. Giblin
# Plot the PDF(SNR) for all redshift bin combinations as a ladder plot

import numpy as np
import pylab as plt
import os
import matplotlib
from matplotlib import rc
from matplotlib import rcParams
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines

# Some lines just to make nice plot fonts
rcParams['ps.useafm'] = False
rcParams['pdf.use14corefonts'] = True
font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 28}
plt.rc('font', **font)
plt.rcParams["mathtext.fontset"] = "cm"

mock_Type = 'cosmoSLICS'
Survey = 'KiDS1000'  #'LSST' #'KiDS1000'
MRres = 'MRres140.64arcs'   # 'MRres60arcs'
Mask = 'Mosaic'        # "Mosaic' 'NoMask'
if Mask == 'Mosaic':
    filetag = 'Mosaic'
else:
    filetag = 'test'

PDForCUM = 'PDF' # Plot PDF or cumul(PDF): "PDF" or "cumPDF"
Plot_Ratio = True # Plot the ratio of measurements to the fiducial cosmol.

Shuffle_Config = 0 # The seed used to shuffle the LOS&R in SLICS

SS=2.816 # 1.408 2.816 5.631 ... #3.11 9.33 18.66
cosmol = []
for i in range(25):
    cosmol.append('%s' %i)
cosmol.append('fid')
#cosmol = 'fid'

if PDForCUM == 'PDF':
    stat_keyword = ''
    plot_keyword = ''
    power = -6. #0.    # 10^scaling applied to statistic in plotting
    
elif PDForCUM == 'cumPDF':
    stat_keyword = 'cum'
    plot_keyword = r'$\Sigma$ '
    power = 0.0    # 10^scaling applied to the statistic in plotting


# Read in PDF's for all zbin combinations
Survey = 'KiDS1000'
if Survey == 'KiDS1000':
    noise = [ 'SN0.27', 'SN0.258', 'SN0.273', 'SN0.254', 'SN0.27']
    noise_KiDS = [ 'SNKiDS', 'SNKiDS', 'SNKiDS', 'SNKiDS', 'SNKiDS']
    # to read in SNKiDS-cosmoSLICS, need to change indir & filename below.
    #noise = [ 'SN0.265' ]
    NLOS = 2170
elif Survey == 'LSST':
    noise = [ 'SN0.28', 'SN0.28', 'SN0.28', 'SN0.28', 'SN0.28']
    #NLOS = 616
ZBcut = ['0.1-0.3', '0.3-0.5', '0.5-0.7', '0.7-0.9', '0.9-1.2']
#ZBcut = ['0.1-1.2']

nkbins = 4 #2000 #4
nzbins_auto = len(noise)
nzbins = nzbins_auto * (nzbins_auto + 1) // 2
PDFs = np.zeros([ len(cosmol), nzbins, nkbins ])

PDFs_data = np.zeros([ nzbins, nkbins ]) # read in n store SLICS data vector
PDFs_KiDS = np.zeros_like( PDFs_data )

# Read in Sys measurements:
Plot_Sys = True
Sys = "IA" # IA / BaryON / dz
if Sys=="IA":
    Sys_vals = [0.0] #[-3.0, -1.0, 0.0, 1.0, 3.0] #np.arange(-6., 7.0) # actual vals to use for colour scale
    Sys_labs = [0.0] #[-3.0, -1.0, 0.0, 1.0, 3.0] #np.arange(-6., 7.0) # what they'll be labelled as in legend
    Sys_files= ['IA0.0'] #['IA-3.0', 'IA-1.0', 'IA0.0', 'IA1.0', 'IA3.0'] #np.arange(-6., 7.0) # differences in filename
    leg_Label = r'$A_{\rm IA}=$'
    c_offset = 0 #2 # helps to get colour of Sys lines visible
    
elif Sys=="dz":
    Sys_vals = [-0.8980,-0.2293,0.,0.2623,0.8254]
    Sys_labs = [-0.90,-0.23,0.,0.26,0.83]
    Sys_files= ['dz5','dz3','dz0','dz4','dz2']  
    leg_Label = r'$\delta_z/\sigma_z=$'
    c_offset = 0 #0.5

elif Sys=="BaryON":
    Sys_vals = [0, 1]
    Sys_labs = ["Baryons", "DM-only"]
    Sys_files= ["BaryON", "BaryOFF"]
    leg_Label = ''
    c_offset = 0 #0.5
    
PDFs_Sys = np.zeros([ len(Sys_vals), nzbins, nkbins ])
PDFs_SLC = np.zeros([ nzbins, nkbins ]) # impact of source-lens clustering given by diff
PDFs_SLC0 = np.zeros([ nzbins, nkbins ])# in IA0 and fid cosmoSLICS (both non-mosaic) - ONLY USED FOR Sys=IA

PDFs_Noise1 = np.zeros([ nzbins, nkbins ]) # Gauss Noise
PDFs_Noise2 = np.zeros([ nzbins, nkbins ]) # KiDS-Noise

# given by the difference of cosmoSLICS-fid and IA0.0
PDFs_err  = np.zeros([ nzbins, nkbins ])

PDFs_fid_mod = np.zeros_like( PDFs_data )  # the fid PDF scaled to include Sys
PDFs_fid_mod2 = np.zeros_like( PDFs_data ) # same, but Sys contribution is additive      

for c in range(len(cosmol)):
    k=0
    for i in range(len(noise)):
        for j in range(i, len(noise)):
            indir = '%s_100Sqdeg_%s_%s_%sGpAM_z%s_ZBcut%s' %(MRres,noise[i],Mask,Survey,Survey,ZBcut[i])
            if j>i:
                indir += '_X_ZBcut%s' %ZBcut[j]
            if mock_Type == 'cosmoSLICS':
                indir += '_Cosmol%s' %cosmol[c]

            filename = indir + '/%s_%s.%sGpAM.LOSAll.SS%s.SNR%sPDF_%sbins.dat' %(noise[i],filetag,Survey,SS, stat_keyword,nkbins)
            x_array, PDFs[c,k,:] = np.loadtxt(filename, usecols=(0,1), unpack=True)

            if c == len(cosmol)-1 and Survey == 'KiDS1000':
                # KiDS-1000:
                filename_KiDS = filename.replace('%s_'%noise[i], 'SNKiDS_').replace('_Cosmol%s' %cosmol[-1], '_CosmolKiDS1000')
                # SLICS OR Trial42 (changed as you prefer)
                filename_data = filename.replace('_Cosmol%s' %cosmol[-1], '').replace('SNKiDS',noise[i])
                #'_CosmolTrial42-NOISE0' )
                filename_data = filename_data.replace('LOSAll', 'LOSAll-Shuffle%s' %Shuffle_Config) # for SLICS only

                # NOISE (Gauss & KiDS)
                filename_Noise1 = filename.replace('%s_'%noise[i], 'NOISE_').replace('_Cosmol%s' %cosmol[-1], '').replace('LOSAll', 'LOSAll-Shuffle0')
                filename_Noise2 = filename_Noise1.replace('NOISE_', 'NOISE-KiDS_')
                
                # Cov:
                filename_cov = filename_data.replace('_%s'%noise[i],'_SNCycle')
                filename_cov = filename_cov.replace('LOSAll', 'NLOS%s' %NLOS)                                                            
                filename_cov = filename_cov.replace('.dat', '.CovMat.npy') 
                
                # read in the Sys vector:
                if Plot_Sys:
                    for s in range(len(Sys_vals)):
                        # Use this to select Sys-contam'd mosaic mocks (perf'd by Make_PDFs_NonMosaic):
                        tmp_bias = '.bias%s-added.dat'%Sys_files[s]
                        filename_sys = filename.replace('.dat', tmp_bias)
                        PDFs_Sys[s,k,:] = np.loadtxt( filename_sys, usecols=(1,), unpack=True )
                    if Sys=="IA" and Sys_labs[s]==0.0:
                        # the IA mocks contain source-lens clustering:
                        filename_slc = filename_sys.replace('MRres140.64arcs_','MRres140.64arcs_SLC1.25_').replace('Mosaic_Ki','NoMask_Ki').replace('Mosaic.K','test.K').replace('bias%s-added.dat'%Sys_files[s],'dat')
                        filename_slc0 = filename_slc.replace('_SLC1.25','')
                        PDFs_SLC[k,:] = np.loadtxt( filename_slc, usecols=(1,), unpack=True )
                        PDFs_SLC0[k,:] = np.loadtxt( filename_slc0, usecols=(1,), unpack=True )
                        
                #print( filename_data )
                PDFs_data[k,:] = np.loadtxt( filename_data, usecols=(1,), unpack=True )
                PDFs_KiDS[k,:] = np.loadtxt( filename_KiDS, usecols=(1,), unpack=True )
                # save an SLC-corrected KiDS measurement?
                if Plot_Sys and Sys=="IA" and Sys_labs[s]==0.0:
                    np.savetxt(filename_KiDS.replace('.dat', '.SLC1.25-corrected.dat'),
                               np.c_[x_array, (PDFs_KiDS-PDFs_SLC+PDFs_SLC0)[k,:]],
                               header='SLC correction applied additively')
                    
                #PDFs_Noise1[k,:] = np.loadtxt( filename_Noise1, usecols=(1,), unpack=True )
                #PDFs_Noise2[k,:] = np.loadtxt( filename_Noise2, usecols=(1,), unpack=True )
                PDFs_err[k,:] = np.sqrt( np.diag( np.load(filename_cov) ) )
            k+=1

#sys.exit()
            
def Ladder_Plot(x_array, y_array,
                y_array_data, y_array_err,
                y_array_data2, 
                y_array_sys, y_array_slc, # probes SLC in IA mocks
                n_row_col, x_label, y_label, x_lims, y_lims, savename):

    fig = plt.figure(figsize = (19,10)) #figsize = (13.5,10)
    gs1 = gridspec.GridSpec(n_row_col, n_row_col)
    d=0     # data number
    
    axes = []
    for j in range(n_row_col):   # scroll across columns 
        for i in range(n_row_col):       # scroll down rows                                                             

            #axes.append( plt.subplot(gs1[p]) )
            ax1 = plt.subplot(gs1[i,j]) #gs1[p])
            if j>i:
                # Dont plot panels above the diagonal                                                                        
                ax1.axis('off')
            else:
                # MAKE cosmoSLICS a FILLED REGION
                upper = np.zeros(len(x_array))
                lower = np.zeros(len(x_array))
                for x in range(len(x_array)):
                    upper[x] = y_array[:,d,x].max()
                    lower[x] = y_array[:,d,x].min()
                plt.fill_between(x_array, lower, upper, color='dimgrey') #, alpha=0.5)
                
                for c in range(y_array.shape[0]):
                    if Plot_Ratio:
                        #if pcosmo[c] in [0.6101, 0.6615, 0.7232, 0.7821, 0.8321, 0.8947]:
                            # plot only the most extreme S8's:
                        ax1.plot( x_array, y_array[c,d,:], color='dimgrey') #color=s_m2.to_rgba(diff_pcosmo[c]) )
                            #continue
                    else:
                        ax1.plot( x_array, y_array[c,d,:], color=s_m.to_rgba(pcosmo[c]) )
                
                handles = []
                # plot Sys vs non-Sys
                if Plot_Sys:
                    for s in range(len(Sys_vals)):
                        ax1.plot( x_array, y_array_sys[s,d,:], color=s_m3.to_rgba(Sys_vals[s]), linewidth=4 )
                        handles.append( mlines.Line2D([], [], color=s_m3.to_rgba(Sys_vals[s]),
                                                      linestyle='-', label=leg_Label+r'%s' %(Sys_labs[s]) ))
                    if Sys=="IA" and Sys_labs[0]==0.0:
                        # plot the SLC measurement
                        s=0 # just the IA0 result
                        ax1.plot( x_array, y_array_slc[d,:], color='pink', linewidth=4, linestyle=':' )
                        handles.append( mlines.Line2D([], [], color='pink', linestyle=':', linewidth=4,label='SLC') )
                    
                        
                #ax1.plot( x_array, y_array[-1,d,:], color='dimgrey', linewidth=2 )
                ax1.errorbar( x_array, y_array_data[d,:], yerr=y_array_err[d,:], color='magenta', linewidth=2 )
                ax1.errorbar( x_array, y_array_data2[d,:], yerr=y_array_err[d,:],
                              color='magenta', linestyle='--', linewidth=2 )
                # if we scale KiDS-PDFs by (Gauss-Noise/KiDS-Noise), does it look more like sims? NO! Sadly.
                #ax1.errorbar( x_array, y_array_data[d,:] * (PDFs_Noise1/PDFs_Noise2)[d,:],
                #              yerr=y_array_err[d,:], color='cyan', linestyle='--',linewidth=2 )
                #ax1.plot( x_array, y_array_data2[d,:], color='cyan', linewidth=2 )
                #ax1.plot( x_array, y_array_data3[d,:], color='dimgrey', linewidth=2 )
                
                #handles.append( mlines.Line2D([], [], color='magenta', linestyle='-', label=r'KiDS-1000' ))
                #handles.append( mlines.Line2D([], [], color='magenta', linestyle='--', label=r'Trial42-NOISE0' ))
                #handles.append( mlines.Line2D([], [], color='cyan', linestyle='-', label=r'$A_{\rm IA}=1.0$' ))
                #handles.append( mlines.Line2D([], [], color='dimgrey', linestyle='-', label=r'SLICS' ))

                ax1.text( 0.60, 0.8, r'%s-%s' %(j+1,i+1),
                                 horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes)
                ax1.set_xlim(x_lims)
                ax1.set_ylim(y_lims)

                # Get rid of x/y ticks completely for subplots not at the edge                                               
                if j==0 and i==2:
                    ax1.set_ylabel(y_label)
                if j>0:
                    ax1.set_yticks([])

                if i==n_row_col-1:
                    ax1.set_xlabel( x_label )
                else:
                    ax1.set_xticks([])
                d+=1

    fig.subplots_adjust(hspace=0, wspace=0)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    if Plot_Ratio:
        fig.colorbar(s_m2, cax=cbar_ax, label=r'$\Delta %s$' %cbar_label)
    else:
        fig.colorbar(s_m, cax=cbar_ax, label=r'$%s$' %cbar_label)

    if Plot_Sys:
        fig.legend(handles=handles, loc='upper center', frameon=False)

    plt.savefig(savename)
    plt.show()
    return

# Read in the S8 values for the cosmoSLICS
ColourBy = 'S8'
if ColourBy == 'S8':
    col=1
    cbar_label = 'S_8'
elif ColourBy == 'h':
    col=2
    cbar_label = 'h'
    
pcosmo = np.loadtxt('/home/bengib/cosmoSLICS/Cosmologies/SLICS_Cosmol_Table_ExtraCol_25plusFid.txt',
                usecols=(1,), unpack=True)
norm = matplotlib.colors.Normalize(vmin=pcosmo.min(), vmax=pcosmo.max())
c_m = matplotlib.cm.jet #cool
s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
s_m.set_array([])

# also make one of these for difference in S_8 relative to fid cosmology:
diff_pcosmo = pcosmo-pcosmo[-1]
norm2 = matplotlib.colors.Normalize(vmin=diff_pcosmo.min(), vmax=diff_pcosmo.max())
s_m2 = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm2)
s_m2.set_array([])

# also make one for colour-coding by Sys:
# setting low lim to IA.min()-2, so no line is WHITE
norm3 = matplotlib.colors.Normalize(vmin=min(Sys_vals)-c_offset, vmax=max(Sys_vals))
c_m3 = matplotlib.cm.cool #Greys
s_m3 = matplotlib.cm.ScalarMappable(cmap=c_m3, norm=norm3)
s_m3.set_array([])

y_limits= None #[-0.05, 1.05]
#sys.exit()

if SS==3.11 or SS==1.408:
    alpha = 0.0 # plotting SNR^alpha x PDF
    x_limit = 3.4 #1.2
    if SS==1.408:
        ss_label = '3.3 arcmin'
    else:
        ss_label = '2.2 arcmin'
    if PDForCUM == 'PDF':
        y_limits=[0.,9.9] #[0.9,2.5] #[0,1.9]
    
elif SS==9.33 or SS==2.816:
    alpha = 0.0
    x_limit = 3.4  #2.3
    ss_label = '6.6 arcmin'
    if PDForCUM == 'PDF':
        y_limits=[0.,9.9]  # for the PDF/1e3

elif SS==18.66 or SS==5.631:
    alpha = 0.0
    x_limit = 3.4  #2.3

    ss_label = '13.2 arcmin'
    if PDForCUM == 'PDF':
        y_limits=[0.,9.9]  #[0,2.9]  # for the PDF/1e3   
    
elif SS == 84.85:
    alpha = 0.0
    x_limit = 2.3
    ss_label = '60.0 arcmin'
    if PDForCUM == 'PDF':
        y_limits=[0,19]


if Plot_Ratio:
    # set limits to appropriate percentage difference:
    y_limits= [-8,8] #[-1.9,1.9] #[-19,19]
    Ratios = PDFs / PDFs[-1]
    # marco suggests dividing by PDF_err (tried it & cant see diffs) BUT ratio looks worse because PDF small in tails
    Ratios_data = PDFs_data / PDFs[-1]
    Ratios_err = PDFs_err / PDFs[-1]
    Ratios_KiDS = PDFs_KiDS / PDFs[-1]
    
    # Get rid of the infinities and nans
    Ratios[np.where(np.isfinite(Ratios) == False)] = 1.
    Ratios_data[np.where(np.isfinite(Ratios_data) == False)] = 1.
    Ratios_err[np.where(np.isfinite(Ratios_err) == False)] = 0.

    if Plot_Sys:
        Ratios_Sys = PDFs_Sys / PDFs[-1]
        Ratios_Sys[np.where(np.isfinite(Ratios_Sys) == False)] = 1.
        if Sys=="IA" and Sys_labs[s]==0.0:
            Ratios_SLC = (PDFs[-1] + PDFs_SLC-PDFs_SLC0) / PDFs[-1]
            Ratios_SLC[np.where(np.isfinite(Ratios_SLC) == False)] = 1.
            Ratios_SLC[np.where(abs(Ratios_SLC) > 1e7)] = 1. # some finite, but massive
            # can we correct the SLC in KiDS?
            PDFs_KiDS_SLC = PDFs_KiDS -PDFs_SLC + PDFs_SLC0 # additive correction
            #PDFs_KiDS_SLC = PDFs_KiDS *PDFs_SLC0/PDFs_SLC # multiplicative correction
            Ratios_KiDS_SLC = PDFs_KiDS_SLC / PDFs[-1]
            
    else:
        Ratios_Sys=0.


    # PLOT THE %-DIFF WITH THE FID COSMOL.
    Ladder_Plot(x_array,
                100*(Ratios-1.), 100*(Ratios_KiDS-1.), 100*Ratios_err,
                100*(Ratios_KiDS_SLC-1.),
                #100*(Ratios_data-1.),
                100*(Ratios_Sys-1.), 100*(Ratios_SLC-1.),
                len(noise), r'SNR',
                r'$\left( {\rm PDF} - {\rm PDF}_{\rm fid} \right)/ {\rm PDF}_{\rm fid}$ '+'[%]'+r' | $\sigma_{\rm s}=$%s'%(ss_label),
                [-1*x_limit, x_limit], y_limits,
                'Figures_4_Paper/Ladder_CS-%sPDF_%s_%s_SS%s_Ratio.png' %(stat_keyword,Survey,Mask,SS) )
else:
    # PLOT THE EXPONENTIAL OF PDFs
    scaling = 10**power
    Ladder_Plot(x_array,
                #np.log10(PDFs), np.log10(PDFs_data), 0.,
                PDFs/PDFs.mean(), PDFs_data/PDFs.mean(), PDFs_err/PDFs.mean(), PDFs_Sys/PDFs.mean(),
                #np.exp(PDFs*scaling), np.exp(PDFs_data*scaling), np.exp(PDFs_data*scaling)*PDFs_err*scaling,
                len(noise), r'SNR',
                r'${\rm PDF} / \langle {\rm PDF} \rangle$ | $\sigma_{\rm s}=$%s' %(ss_label),
                #plot_keyword + r'$\exp\left[ {\rm PDF} \times 10^{%s} \right]$ | $\sigma_{\rm s}=$%s'%(power, ss_label),
                [-1*x_limit, x_limit], y_limits,
                'Figures_4_Paper/Ladder_CS-%sPDF_%s_%s_SS%s.png' %(stat_keyword,Survey,Mask,SS) )


