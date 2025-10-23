# 03/12/2024, B. Giblin
# This reads in the sys. biases for dz or baryons, made by Make_SNR_SNR-PDF.py
# and computes a linear relation as a func of dz or baryons ON/OFF, for each
# zbin and theta-bin separately. This linear relation will be used to model
# the sys impact in the l'hood sampling code.
# -----
# Note this can also do it for IAs, but it's not necessary, since we have
# more A_IA values than dz/baryon data points, and so we've built a good
# GP emulator for this systematic.

import numpy as np
import pylab as plt
import glob
import os
import matplotlib
from matplotlib import rc
from matplotlib import rcParams
import matplotlib.gridspec as gridspec

# Some lines just to make nice plot fonts
rcParams['ps.useafm'] = True
rcParams['pdf.use14corefonts'] = True
font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 20}
plt.rc('font', **font)
plt.rcParams["mathtext.fontset"] = "cm"

#Plot_Ratio = True      # plot the ratio of the measurements to the fid. cosmology.
Sys = "dz"
if Sys=="IA":
    Sys_vals = np.arange(-6., 7.0) # actual vals to use for colour scale
    Sys_labs = Sys_vals # what they'll be labelled as in legend
    Sys_files= ['IA-6.0','IA-5.0','IA-4.0','IA-3.0','IA-2.0','IA-1.0','IA0.0',
                'IA1.0','IA2.0','IA3.0','IA4.0','IA5.0','IA6.0'] #np.arange(-6., 7.0) # differences in filename
    leg_Label = r'$A_{\rm IA}$'
    c_offset = 2 #2 # helps to get colour of Sys lines visible     
    
    NLOS_Sys = 50   
    Sys_Tag0 = 'IA0.0'
    Mask_Sys = 'NoMask'                                                        
    name_Sys = 'test'

elif Sys=="dz":
    Sys_vals = np.array([-0.8980,-0.2293,0.2623,0.8254])
    Sys_labs = [-0.90,-0.23,0.26,0.83]
    Sys_files= ['dz5','dz3','dz4','dz2']
    leg_Label = r'$\delta_z$'
    c_offset = 0.5 #0.5 

    NLOS_Sys = 18
    Sys_Tag0 = 'dz5' # use dz5 as benchmark in dz case.
    Mask_Sys = 'Mosaic'
    name_Sys = 'Mosaic'

elif Sys=="BaryON":
    Sys_vals = np.array([0, 1])
    Sys_labs = ["DM-only", "Baryons",]
    Sys_files= ["BaryOFF", "BaryON"]
    leg_Label = 'Baryons'
    c_offset = 0.5
                    
    NLOS_Sys = 180
    Sys_Tag0 = 'BaryOFF'
    Mask_Sys = 'Mosaic'
    name_Sys = 'Mosaic'

cosmol = ['fid']

# Read in PDF's for all zbin combinations
MRres = '140.64arcs' # '60arcs'
Mask = 'Mosaic'
SS = 2.816           # smoothing scale used in clipped measurements

if Mask == 'Mosaic':
    filetag = 'Mosaic'
    NLOS_CS = 900
else:
    filetag = 'test'
    NLOS_CS = 50

Survey = 'KiDS1000'
if Survey == 'KiDS1000':
    noise = ['SN0.27', 'SN0.258', 'SN0.273', 'SN0.254', 'SN0.27'] #['SN0.265']
    NLOS = 2170
elif Survey == 'LSST':
    noise = [ 'SN0.28', 'SN0.28', 'SN0.28', 'SN0.28', 'SN0.28']
    NLOS = 616
ZBcut = ['0.1-0.3', '0.3-0.5', '0.5-0.7', '0.7-0.9', '0.9-1.2'] #['0.1-1.2']
    
nthbins = 4
nzbins = 15

pdf_CS_Sys_bias = np.zeros([ len(Sys_vals), nzbins, nthbins ]) # the pure additive bias (bias_x-bias_0)
'''
pdf_err  = np.zeros([ nzbins, nthbins ])                    # SLICS error
'''

# define function for linear fit (just grad; defined so bias=0 at x[0])
from scipy.optimize import curve_fit
# linear model used to fit sys-bias per bin (e.g. per theta for 2pcf) as a func of Sys_Nodes x (e.g. dz).
# Because curve_fit doesnt play nice with fixed params (x0), need an outer & an inner loop to provide this. 
def linear(x0):
    def inner_linear(x, m):
        return m*(x-x0)
    return inner_linear

# store the gradient (hyperparam) & its error for each z&theta bin
HPs = np.zeros([nzbins, nthbins ]) # unclipped
HPs_err = np.zeros([nzbins, nthbins ])

for c in range(len(cosmol)):
    k=0
    print("On cosmology ", cosmol[c])
    for i in range(len(noise)):
        for j in range(i, len(noise)):
            indir_prefix = 'MRres%s_%s_100Sqdeg_%s_%s_%sGpAM_z%s_' %(MRres,Sys_Tag0,noise[i],Mask_Sys,Survey,Survey)
            if j != i:
                ZBlabel = 'ZBcut%s_X_ZBcut%s' %(ZBcut[i],ZBcut[j])
            else:
                ZBlabel = 'ZBcut%s' %(ZBcut[i])
            indir_Sys = '%s%s_Cosmol%s' %(indir_prefix,ZBlabel,cosmol[c])
            filename_Sys = indir_Sys + '/%s_%s.%sGpAM.LOSAll.SS%s.SNRPDF_%sbins.bias%s-additive.dat' %(noise[i],name_Sys,
                                                                                                       Survey,SS,nthbins,Sys_Tag0)
            if c == len(cosmol)-1:
                # ---------------------------------- Sys READING ------------------------------------------
                for s in range(len(Sys_vals)):
                    # unclipped:
                    filename_Sys2 = filename_Sys.replace(Sys_Tag0, Sys_files[s])

                    if Sys == "IA":
                        tmp_bias = 'IAbias%s'%Sys_vals[s]
                    elif Sys == "dz" or Sys == "BaryON":
                        tmp_bias = 'bias%s'%Sys_files[s]
                    
                    # now load the pure bias:
                    SNR, pdf_CS_Sys_bias[s,k,:] = np.loadtxt( filename_Sys2, usecols=(0,1),unpack=True)
                # fit a linear relation for each z&SNR bin 
                for t in range(nthbins):
                    HPs[k,t], cov = curve_fit(linear(x0=Sys_vals[0]), Sys_vals, SNR[t]*pdf_CS_Sys_bias[:,k,t])

                # In Make_SNR-PDF, the dz bias was computed relative to dz=-0.89; so bias=0 for this point.
                # We need to recalibrate this, so bias=0 for dz=0, as it is in the mocks & data.
                # We can't use the simple SLICS mocks for a dz=0 pred; tried this, it's outside of the dz2-5 results (for some reason).
                # So we now use the fits to compute the bias @ dz=0, and subtract this off, saving new
                # pure biases to be used in the MCMC, which will be centred on 0.
                if Sys=="dz":
                    b0 = np.zeros(nthbins)
                    for t in range(nthbins):
                        b0[t] = linear(x0=Sys_vals[0])(0., HPs[k,t]) / SNR[t]
                        # bias at dz=0 given by old fit^ (units of PDF)
                    
                    for s in range(len(Sys_files)):
                        # subtract bias(dz=0) from previous pure biases:
                        new_bias = pdf_CS_Sys_bias[s,k,:] - b0
                        # save the new bias:                                                                                      
                        new_fname_bias = filename_Sys.replace(Sys_Tag0, Sys_files[s]).replace('additive.dat',
                                                                                              'additive.0centred.dat')
                        print(new_fname_bias)
                        print("--------")
                        np.savetxt(new_fname_bias, np.c_[SNR, new_bias], header='# SNR, dz-bias[rel. to dz=0]')
                        # replace the bias:
                        pdf_CS_Sys_bias[s,k,:] = new_bias
                
            k+=1

# Load the SLICS err corresponding to a complete (shuffled) survey footprint
# this will be used as the error on the Sys measurement:
'''
indir_S = 'MRres%s_100Sqdeg_SNCycle_%s_%sGpAM_z%s_ZBcut0.1-0.3_X_ZBcut0.1-0.3/ThBins%s/NLOS%s/' %(MRres,Mask,Survey,Survey,
                                                                                                  nthbins,NLOS)
# unclipped:
filename_cov = '%s_%s.%sGpAM.NLOS%s-Shuffle0.ORIG.UCxUC%s.CovMat.zbins11_12_13_14_15_22_23_24_25_33_34_35_44_45_55.npy' %(noise[i],filetag,Survey,NLOS,pm)
filename_cov = indir_S + filename_cov
xi_err = np.sqrt( np.diag( np.load(filename_cov) ) ) # stdev for one K1000 survey                               
xi_err = np.reshape(xi_err, (nzbins,nthbins))        # reshape
xi_err = xi_err / np.sqrt(NLOS_Sys/18.)              # make it error on the mean of 1 survey
'''

def Ladder_Plot(x_array,
                y_array, y_array_err, hps, # hyperparams (grad) for linear fit
                n_row_col, x_label, y_label, x_lims, y_lims, savename):

    fig = plt.figure(figsize = (13.5,10)) #figsize = (20,14)
    gs1 = gridspec.GridSpec(n_row_col, n_row_col)
    d=0     # data number
    
    axes = []
    for j in range(n_row_col):   # scroll across columns 
        for i in range(n_row_col):       # scroll down rows                                                                   
            ax1 = plt.subplot(gs1[i,j])
            if j>i:
                # Dont plot panels above the diagonal                                                                        
                ax1.axis('off')
            else:
                for t in range(nthbins):
                    #ax1.plot( x_array, y_array[:,d,t], color=s_m.to_rgba(np.log10(theta[t])) )
                    ax1.plot( x_array, linear(x0=0.)(x_array, hps[d,t]), linestyle='-', color=s_m.to_rgba(SNR[t]) )
                    ax1.errorbar( x_array, y_array[:,d,t], yerr=0., #y_array_err[d,t],
                                  fmt='o', ecolor=s_m.to_rgba(SNR[t]),
                                  markerfacecolor=s_m.to_rgba(SNR[t]),
                                  markeredgecolor=s_m.to_rgba(SNR[t]) )
                    
                ax1.text( 0.75, 0.8, r'%s-%s' %(j+1,i+1),
                                 horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes)
                #ax1.set_yscale('log')
                #ax1.set_xscale('log')
                ax1.set_xlim(x_lims)
                if len(y_lims) > 2:
                    ax1.set_ylim(y_lims[i])
                else:
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
    fig.colorbar(s_m, cax=cbar_ax, label=r'SNR')
    #plt.savefig(savename)
    plt.show()
    return

# also make a colour scale for theta
c_m = matplotlib.cm.jet
norm = matplotlib.colors.Normalize(vmin=SNR.min(), vmax=SNR.max()) 
s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
s_m.set_array([])

if SS==3.11:
    ss_label = '3.3 arcmin'
elif SS==2.816:
    ss_label = '6.6 arcmin'
elif SS==5.631:
    ss_label = '13.2 arcmin'

y_limits = [-9., 11.]
    
Ladder_Plot(Sys_vals,
            SNR*pdf_CS_Sys_bias, 0.,
            HPs, len(noise), leg_Label,
            r'$\Delta {\rm PDF}$ | $\sigma_{\rm s}=$%s' %(ss_label),
            None, y_limits, None)


