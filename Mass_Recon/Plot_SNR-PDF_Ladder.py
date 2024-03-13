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
Plot_Ratio = False # Plot the ratio of measurements to the fiducial cosmol.

A_IA = 0.0  # Intrinsic alignment amplitude
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
    #noise = [ 'SN0.27', 'SN0.258', 'SN0.273', 'SN0.254', 'SN0.27']
    noise = [ 'SN0.265' ]
    #NLOS = 3906 #715
elif Survey == 'LSST':
    noise = [ 'SN0.28', 'SN0.28', 'SN0.28', 'SN0.28', 'SN0.28']
    #NLOS = 616
#ZBcut = ['0.1-0.3', '0.3-0.5', '0.5-0.7', '0.7-0.9', '0.9-1.2']
ZBcut = ['0.1-1.2']

nkbins = 2000 #2000 #4
nzbins_auto = len(noise)
nzbins = nzbins_auto * (nzbins_auto + 1) // 2
PDFs = np.zeros([ len(cosmol), nzbins, nkbins ])

PDFs_data = np.zeros([ nzbins, nkbins ])
PDFs_data2 = np.zeros([ nzbins, nkbins ])
PDFs_data3 = np.zeros([ nzbins, nkbins ])
PDFs_err  = np.zeros([ nzbins, nkbins ])

PDFs_fid_mod = np.zeros_like( PDFs_data )  # the fid PDF scaled to include IA
PDFs_fid_mod2 = np.zeros_like( PDFs_data ) # same, but IA contribution is additive      

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
                filename_data = filename.replace('_Cosmol%s' %cosmol[-1], '')
                filename_data = filename_data.replace('LOSAll', 'LOSAll-Shuffle%s' %Shuffle_Config)
                # read in the IA vector:
                #filename_data = indir.replace('%s_'%MRres, '%s_IA%s_'%(MRres,A_IA)) + '/%s_%s.%sGpAM.LOSAll.SS%s.SNR%sPDF_%sbins.dat' %(noise[i],filetag,Survey,SS, stat_keyword,nkbins)
                #print( filename_data )
                PDFs_data[k,:] = np.loadtxt( filename_data, usecols=(1,), unpack=True )

                #filename_data2 = filename_data.replace('_IA%s_'%A_IA, '_IA1.0_')
                #PDFs_data2[k,:] = np.loadtxt( filename_data2, usecols=(1,), unpack=True )
                
                #filename_cov = filename.replace('SNCycle', noise[i]) # DONT NEED THIS NOW NOT DOING SNCycle
                #filename_cov = filename.replace('_Cosmolfid', '')
                # read in SLICS data vector
                #filename_data3 = filename_cov
                #PDFs_data3[k,:] = np.loadtxt( filename_data3, usecols=(1,), unpack=True )
                
                #filename_cov = filename_cov.replace('LOSAll', 'NLOS%s' %NLOS)
                #filename_cov = filename_cov.replace('.dat', '.CovMat.npy')
                #print( filename_cov )
                PDFs_err[k,:] = np.zeros(nkbins) #np.sqrt( np.diag( np.load(filename_cov) ) / 10. ) # scaled from 100 --> 1000 deg^2

                # manually add the IA contribution to the fid cosmology and save it:
                #PDFs_fid_mod[k,:] = PDFs[-1,k,:] * ( PDFs_data2[k,:]/PDFs_data[k,:] )
                #PDFs_fid_mod2[k,:] = PDFs[-1,k,:] + ( PDFs_data2[k,:]-PDFs_data[k,:] )
                # save the modified data
                #filename_mod = filename.replace('%sbins.dat' %nkbins, '%sbins.IA1.0-scaled.dat' %nkbins)
                #filename_mod2 = filename.replace('%sbins.dat' %nkbins, '%sbins.IA1.0-added.dat' %nkbins)
                #print( filename_mod )
                #print(filename_mod2)
                #print("----------------")
                #np.savetxt( filename_mod, np.c_[x_array,PDFs_fid_mod[k,:]] )
                #np.savetxt( filename_mod2, np.c_[x_array,PDFs_fid_mod2[k,:]] )
                
                
            k+=1

#sys.exit()
            
def Ladder_Plot(x_array, y_array,
                y_array_data, y_array_err,
                #y_array_data2, y_array_data3,
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
                for c in range(y_array.shape[0]):
                    if Plot_Ratio:
                        #if pcosmo[c] in [0.6101, 0.6615, 0.7232, 0.7821, 0.8321, 0.8947]:
                            # plot only the most extreme S8's:
                        ax1.plot( x_array, y_array[c,d,:], color=s_m2.to_rgba(diff_pcosmo[c]) )
                            #continue
                    else:
                        ax1.plot( x_array, y_array[c,d,:], color=s_m.to_rgba(pcosmo[c]) )

                # plot IA vs non-IA
                #ax1.plot( x_array, y_array[-1,d,:], color='dimgrey', linewidth=2 )
                #ax1.plot( x_array, y_array_data[d,:], color='magenta', linewidth=2 )
                #ax1.plot( x_array, y_array_data2[d,:], color='cyan', linewidth=2 )
                #ax1.plot( x_array, y_array_data3[d,:], color='dimgrey', linewidth=2 )
                
                handles = []
                #handles.append( mlines.Line2D([], [], color='magenta', linestyle='-', label=r'$A_{\rm IA}=0.0$' ))
                #handles.append( mlines.Line2D([], [], color='cyan', linestyle='-', label=r'$A_{\rm IA}=1.0$' ))
                handles.append( mlines.Line2D([], [], color='dimgrey', linestyle='-', label=r'SLICS' ))
                
                #ax1.errorbar( x_array, y_array_data[d,:], yerr=y_array_err[d,:]/1e4, color='magenta', linewidth=2 )

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

    # if plotting the IA
    #fig.legend(handles=handles, loc='upper center', frameon=False)

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

y_limits= None #[-0.05, 1.05]

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
    y_limits=[-19,19]
    Ratios = np.zeros_like(PDFs)
    Ratios_data = np.zeros_like(PDFs_data)
    Ratios_data2 = np.zeros_like(PDFs_data)
    Ratios_data3 = np.zeros_like(PDFs_data)
    Ratios_mod   = np.zeros_like(PDFs_data)
    Ratios_err  = np.zeros_like(PDFs_err)

    k = 0
    for i in range(len(noise)):
        for j in range(i, len(noise)):
            for c in range(PDFs.shape[0]):
                Ratios[c,k,:] = PDFs[c,k,:] / PDFs[-1,k,:]
            Ratios_data[k,:] = PDFs_data[k,:] / PDFs[-1,k,:]
            #Ratios_data2[k,:] = PDFs_data2[k,:] / PDFs[-1,k,:]
            #Ratios_data3[k,:] = PDFs_data3[k,:] / PDFs[-1,k,:]
            #Ratios_mod[k,:]   = PDFs_fid_mod2[k,:] / PDFs[-1,k,:]
            Ratios_err[k,:] = PDFs_err[k,:] / PDFs[-1,k,:]
            k += 1

    # Get rid of the infinities and nans
    Ratios[np.where(np.isfinite(Ratios) == False)] = 1.
    Ratios_data[np.where(np.isfinite(Ratios_data) == False)] = 1.
    #Ratios_data2[np.where(np.isfinite(Ratios_data2) == False)] = 1.
    #Ratios_data3[np.where(np.isfinite(Ratios_data3) == False)] = 1.
    #Ratios_mod[np.where(np.isfinite(Ratios_mod) == False)] = 1.
    Ratios_err[np.where(np.isfinite(Ratios_err) == False)] = 0.


    # PLOT THE %-DIFF WITH THE FID COSMOL.
    Ladder_Plot(x_array,
                100*(Ratios-1.), 100*(Ratios_data-1.), 100*Ratios_err,
                #100*(Ratios_data2-1.), 100*(Ratios_data3-1.),
                len(noise), r'SNR',
                r'$\left( {\rm PDF} - {\rm PDF}_{\rm fid} \right)/ {\rm PDF}_{\rm fid}$ '+'[%]'+r' | $\sigma_{\rm s}=$%s'%(ss_label),
                [-1*x_limit, x_limit], y_limits,
                'Figures_4_Paper/Ladder_CS-%sPDF_%s_%s_SS%s_Ratio.png' %(stat_keyword,Survey,Mask,SS) )
else:
    # PLOT THE EXPONENTIAL OF PDFs
    scaling = 10**power
    Ladder_Plot(x_array,
                #np.log10(PDFs), np.log10(PDFs_data), 0.,
                PDFs/PDFs.mean(), PDFs_data/PDFs.mean(), PDFs_err/PDFs.mean(),
                #np.exp(PDFs*scaling), np.exp(PDFs_data*scaling), np.exp(PDFs_data*scaling)*PDFs_err*scaling,
                len(noise), r'SNR',
                r'${\rm PDF} / \langle {\rm PDF} \rangle$ | $\sigma_{\rm s}=$%s' %(ss_label),
                #plot_keyword + r'$\exp\left[ {\rm PDF} \times 10^{%s} \right]$ | $\sigma_{\rm s}=$%s'%(power, ss_label),
                [-1*x_limit, x_limit], y_limits,
                'Figures_4_Paper/Ladder_CS-%sPDF_%s_%s_SS%s_Exp.eps' %(stat_keyword,Survey,Mask,SS) )


