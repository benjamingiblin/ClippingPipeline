# 28/03/2025: Does it matter that the shape noise we generated was Gaussian
# rather than random rotations of the KiDS ellipticities? Let's see.

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
# Some lines just to make nice plot fonts
from matplotlib import rcParams
rcParams['ps.useafm'] = False
rcParams['pdf.use14corefonts'] = True
font = {'family' : 'serif',
	'weight' : 'normal',
        'size'   : 14}
plt.rc('font', **font)
plt.rcParams["mathtext.fontset"] = "cm"

Calc_Kappa = False # if true, read in all kappa maps, make a kappa PDF
                   # if False, read in the pre-made SNR map
Calc_Combine_zbins = True

noise = ['SN0.27', 'SN0.258', 'SN0.273', 'SN0.254', 'SN0.27']
ZBcut = ['0.1-0.3', '0.3-0.5', '0.5-0.7', '0.7-0.9', '0.9-1.2']
SS = 2.816

Nr = 18 # regions
Nlos = 50

nkbins = 4
if Calc_Combine_zbins:
    nzbins = 15
    j_upper = len(noise)
else:
    nzbins = len(noise)
    # j_upper is updated in loop if just doing auto bins

if Calc_Kappa:
    x_edges = np.linspace(-0.04, 0.04, nkbins+1) # this is kappa
    x = x_edges[:-1] + (x_edges[1]-x_edges[0])/2. 
    # gaussian SN
    pdf_g = np.zeros([ nzbins, Nr*Nlos, nkbins ])
    # KiDS SN
    pdf_k = np.zeros([ nzbins, Nr*Nlos, nkbins ])

# gaussian
pdf_g_avg = np.zeros([ nzbins, nkbins ])
pdf_g_std = np.zeros_like(pdf_g_avg)
# KiDS
pdf_k_avg = np.zeros([ nzbins, nkbins ])
pdf_k_std = np.zeros_like(pdf_g_avg)
# the actual K1000 data PDFs
pdf_k1000_avg = np.zeros([ nzbins, nkbins ])

# read in the SLICS cov which measures the cosmic variance:
std_slics = np.zeros([ nzbins, nkbins ])
# read in all the cosmoSLICS (Gauss SN & KiDS SN):
ncos = 25
pdf_g_cs = np.zeros([ ncos, nzbins, nkbins ])
pdf_k_cs = np.zeros_like( pdf_g_cs )

zcount=0
for i in range( len(noise) ):
    if not Calc_Combine_zbins: j_upper=i # only access auto bins
    for j in range(i, j_upper):
        print(" ---- zbin %s-%s --- " %(i+1,j+1) )
        
        dir_g = 'MRres140.64arcs_100Sqdeg_%s_Mosaic_KiDS1000GpAM_zKiDS1000_ZBcut%s_Cosmolfid/' %(noise[i],ZBcut[i])
        if j!=i: dir_g = dir_g.replace('ZBcut%s'%ZBcut[i], 'ZBcut%s_X_ZBcut%s'%(ZBcut[i],ZBcut[j]))

        if Calc_Kappa:
            lcount=0
            for los in range(1,Nlos+1):
                for R in range(1,Nr+1):
                    infile_g = dir_g + '%s_Mosaic.KiDS1000GpAM.LOS%sR%s.SS%s.Ekappa.npy' %(noise[i],los,R,SS)
                    infile_k = infile_g.replace(noise[i],'SNKiDS')
                    kappa_g = np.load(infile_g)
                    kappa_k = np.load(infile_k)
                
                    pdf_g[zcount,lcount,:],_ = np.histogram(kappa_g[kappa_g != 0], bins=x_edges)
                    pdf_k[zcount,lcount,:],_ = np.histogram(kappa_k[kappa_k != 0], bins=x_edges)
                    count+=1
            pdf_g_avg[zcount,:] = np.mean(pdf_g[zcount,:,:], axis=0)
            pdf_k_avg[zcount,:] = np.mean(pdf_k[zcount,:,:], axis=0)
            pdf_g_std[zcount,:] = np.std(pdf_g[zcount,:,:], axis=0)
            pdf_k_std[zcount,:] = np.std(pdf_k[zcount,:,:], axis=0)

        else:
            # read pre-made SNR PDFs (from Make_SNR-PDF_...py)
            infile_g = dir_g + '%s_Mosaic.KiDS1000GpAM.LOSAll.SS%s.SNRPDF_%sbins.dat' %(noise[i],SS,nkbins)
            x, pdf_g_avg[zcount] = np.loadtxt(infile_g, usecols=(0,1), unpack=True)
            infile_k = infile_g.replace(noise[i],'SNKiDS')
            _, pdf_k_avg[zcount] = np.loadtxt(infile_k, usecols=(0,1), unpack=True)
            # read K1000 PDFs:
            infile_k1000 = infile_k.replace('_Cosmolfid','_CosmolKiDS1000')
            _, pdf_k1000_avg[zcount] = np.loadtxt(infile_k1000, usecols=(0,1), unpack=True)
            
            # get cov (first std measured across cosmoSLICS)...:
            cfile_g =  dir_g + '%s_Mosaic.KiDS1000GpAM.NLOS%s.SS%s.SNRPDF_%sbins.CovMat.npy' %(noise[i],Nlos,SS,nkbins)
            pdf_g_std[zcount] = np.sqrt( np.diag( np.load(cfile_g) ) )
            cfile_k = cfile_g.replace(noise[i],'SNKiDS')
            pdf_k_std[zcount] = np.sqrt( np.diag( np.load(cfile_k) ) )
            # ...next read the SLICS cov:
            covdir = dir_g.replace('_Cosmolfid', '')
            covfile = covdir + '%s_Mosaic.KiDS1000GpAM.NLOS217-Shuffle0.SS%s.SNRPDF_%sbins.CovMat.npy' %(noise[i],SS,nkbins)
            print(covfile)
            std_slics[zcount] = np.sqrt( np.diag( np.load(covfile) ) )

            # get all of cosmoSLICS:
            for c in range(ncos):
                pdf_g_cs[c,zcount,:] = np.loadtxt( infile_g.replace('fid',str(c)), usecols=(1,), unpack=True )
                pdf_k_cs[c,zcount,:] = np.loadtxt( infile_k.replace('fid',str(c)), usecols=(1,), unpack=True )
                
        zcount+=1 # increment zbin



# PLOTTING STUFF
if Calc_Kappa:
    # then we've averaged all LOS and R's, so err on mean is:
    mean_err_g = pdf_g_std/np.sqrt(Nlos*Nr)
    mean_err_k = pdf_k_std/np.sqrt(Nlos*Nr)
    xlabel = 'kappa'
else:
    # the std corresponds to a whole survey of 18R, so the avg is only over LOS:
    mean_err_g = pdf_g_std/np.sqrt(Nlos)
    mean_err_k = pdf_k_std/np.sqrt(Nlos)
    xlabel = 'SNR'

# get the cosmoSLICS S8:
S8 = np.loadtxt('/home/bengib/cosmoSLICS/Cosmologies/SLICS_Cosmol_Table_ExtraCol_25plusFid.txt',
                usecols=(1,), unpack=True)
cbar_label = r'$S_8$'
norm = matplotlib.colors.Normalize(vmin=S8.min(), vmax=S8.max())
c_m = matplotlib.cm.jet
s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
s_m.set_array([])
    

def Plot_Auto_Bins():
    fig = plt.figure(figsize = (15,10))
    gs1 = gridspec.GridSpec(1,len(noise))
    for i in range(len(noise)):
        ax1 = plt.subplot(gs1[0,i])
        if Plot_Ratio:
            # plot 50% of cosmic variance
            cosvar = (std_slics/pdf_g_avg)[i]
            ax1.fill_between(x, 1.-fvar*cosvar, 1.+fvar*cosvar, alpha=0.5, color='dimgrey')
            #for c in range(ncos):
            #    ax1.plot( x, pdf_g_cs[c,i]/pdf_g_avg[i], color=s_m.to_rgba(S8[c]) )
                
            ax1.errorbar( x, ratio[i], yerr=ratio_err[i], color='black', label='SN Gauss V KiDS')
            ax1.plot(x, np.ones_like(x), 'k:')
            ylabel=r'$PDF^{KiDS} / PDF^{Gauss}$'
        else:
            ax1.fill_between(x, pdf_g_avg[i]-fvar*std_slics[i], pdf_g_avg[i]+fvar*std_slics[i], alpha=0.5, color='dimgrey')
            ax1.errorbar( x, pdf_g_avg[i], yerr=pdf_g_std[i]/np.sqrt(Nlos*Nr), color='red', label='Gauss SN' )
            ax1.errorbar( x, pdf_k_avg[i], yerr=pdf_k_std[i]/np.sqrt(Nlos*Nr), color='blue', label='KiDS SN' )
            ylabel='PDF'
            #ax1.set_yscale('log')
            
        ax1.set_xlabel(xlabel)
        ax1.text( 0.75, 0.6, r'%s-%s' %(i+1,i+1),
                                     horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes)
        ax1.set_ylim(y_lims)
        if i==0:
            ax1.set_ylabel(ylabel)
            ax1.legend(loc='upper right', frameon=False)
        else:
            ax1.set_yticks([])
    fig.subplots_adjust(hspace=0, wspace=0)
    #fig.subplots_adjust(right=0.8)
    #cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    #fig.colorbar(s_m, cax=cbar_ax, label=r'$S_8$')
    plt.show()
    return





def Ladder_Plot(x_array, y_array,
                y_array_data, y_array_err, y_array_data2, cosvar,
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

                ax1.fill_between(x, 1.-fvar*cosvar[d], 1.+fvar*cosvar[d], alpha=0.5, color='dimgrey')
                for c in range(y_array.shape[0]):
                    ax1.plot( x_array, y_array[c,d,:], color=s_m.to_rgba(S8[c]) )
                # plot the SN comparison:
                if Plot_Ratio:
                    ax1.plot(x, np.ones_like(x), 'k:') # horizontal line

                label1 = 'SN Gauss V KiDS'
                label2 = 'K1000'
                #label1 = 'SN$=$Gauss'
                #label2 = 'SN$=$K1000'
                ax1.errorbar( x, y_array_data[d], yerr=y_array_err[d], color='black', linewidth=2, label=label1)
                #ax1.errorbar( x, y_array_data2[d], yerr=y_array_err[d], color='grey', linestyle='--', linewidth=2, label=label2)
                
                ax1.text( 0.60, 0.8, r'%s-%s' %(j+1,i+1),
                                 horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes)
                ax1.set_xlim(x_lims)
                ax1.set_ylim(y_lims)

                #if j==0 and i==0:
                #    ax1.legend(loc='upper right', frameon=False)
                    
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
    fig.colorbar(s_m, cax=cbar_ax, label=cbar_label)

    #fig.legend(loc='upper center', frameon=False)
    #plt.savefig(savename)
    plt.show()
    return

fvar = 1.0 # fraction of cosmic variance to plot above/below.
Plot_Ratio = True

# take ratio; decide if Gauss SN or KiDS SN is used on the numerator:
baseline = 'kids' #'gauss'
if baseline == 'gauss':
    num = pdf_k_avg
    dnum = pdf_g_avg
    numerr = mean_err_k
    dnumerr = mean_err_g
    y_label_ratio = r'$PDF^{\rm KiDS} / PDF^{\rm Gauss}$'
    ratio_array = pdf_g_cs/dnum
else:
    num = pdf_g_avg
    dnum = pdf_k_avg
    numerr = mean_err_g
    dnumerr = mean_err_k
    y_label_ratio = r'$PDF^{\rm Gauss} / PDF^{\rm KiDS}$'
    ratio_array = pdf_k_cs/dnum
    
ratio = num / dnum
ratio_err = np.sqrt( (num/dnum**2)**2. *dnumerr**2. + (1/dnum)**2 * numerr**2 )
ratio_k1000 = pdf_k1000_avg / dnum


if Plot_Ratio:
    y_array = ratio_array
    y_array_data = ratio
    y_array_data2 = ratio_k1000
    y_array_err = ratio_err
    y_array_data2 = ratio_k1000
    y_lims=[0.96,1.04]
    cosvar = (std_slics/dnum)
    y_label = y_label_ratio
    
else:
    y_array = pdf_g_cs
    y_array_data = pdf_k_avg
    y_array_data2 = pdf_k1000_avg
    y_array_err = pdf_k_std/np.sqrt(Nlos*Nr)
    y_lims=[0, 1.1*pdf_g_avg.max()]
    cosvar = std_slics
    y_label = r'${\rm PDF(SNR)}$'
    
if Calc_Combine_zbins:
    Ladder_Plot(x, y_array,
                y_array_data, y_array_err, y_array_data2, cosvar,
                len(noise), r'SNR', y_label,
                None, y_lims, None)

    # see if K1000 is closer to SNGauss or SNKiDS (it's proper close to both).
    #Ladder_Plot(x, y_array,
    #        pdf_k1000_avg/pdf_g_avg, np.zeros_like(y_array_err), pdf_k1000_avg/pdf_k_avg, cosvar,
    #        len(noise), r'SNR', r'$PDF^{K1000} / PDF^{SN}$',
    #       None, [0.85,1.15], None)
