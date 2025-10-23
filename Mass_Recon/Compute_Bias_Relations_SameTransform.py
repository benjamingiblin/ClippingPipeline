# 03/12/2024, B. Giblin

# GIBLIN; THIS IS A VARIATION ON COMPUTE_BIAS WHERE THE AUTO & X BINS HAVE DIFFERENT TRANSFORMATIONS
# AUTO HAS SNR*BIAS TRANSFORM, CROSS HAS NORM. TO MEAN0 & STD1

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
import sys
import matplotlib
from matplotlib import rc
from matplotlib import rcParams
import matplotlib.gridspec as gridspec
from sklearn import preprocessing # used to normalise data pre-GP

# Get 2D interp func - used for dz sys
sys.path.insert(0, '/home/bengib/GPR_Emulator')
from GPR_Classes import GPR_Emu


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
    Sys_files = ['dz5','dz3','dz4','dz2']
    leg_Label = r'$\delta_z$'
    c_offset = 0.5 #0.5 
    
    Nsys = len(Sys_files)
    # We need a 2nd set of arrays for dz; one to use with X-zbins
    # where dz's can differ
    Sys_vals_all = []
    Sys_files_all = []
    for s1 in range(Nsys):
    	for s2 in range(Nsys):
    		Sys_vals_all.append([Sys_vals[s1],Sys_vals[s2]])
    		if s1==s2:
    			Sys_files_all.append(Sys_files[s1])
    		else:
    			Sys_files_all.append('%s_X_%s' %(Sys_files[s1],Sys_files[s2]))
    Sys_vals_all = np.array(Sys_vals_all)

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

pdf_Sys_bias = []

# define function for linear fit (just grad; defined so bias=0 at x[0])
from scipy.optimize import curve_fit
# linear model used to fit sys-bias per bin (e.g. per theta for 2pcf) as a func of Sys_Nodes x (e.g. dz).
# Because curve_fit doesnt play nice with fixed params (x0), need an outer & an inner loop to provide this. 
def linear(x0):
    def inner_linear(x, m):
        return m*(x-x0)
    return inner_linear

# store the gradient (hyperparam) & its error for each z&theta bin
#HPs = np.zeros([nzbins, nthbins ]) 
#HPs_err = np.zeros([nzbins, nthbins ])
HPs = [] # Only used for dz sys; same as HPs for auto-bins (lin fit), but differs for X-bins (where a GP is used)

k=0
for i in range(len(noise)):
    for j in range(i, len(noise)):
        # the num of sys-params differs for dz-cross bins, hence:
        if Sys=="dz" and i!=j:
            pdf_Sys_bias_perz = np.zeros([len(Sys_vals_all),nthbins])
            sys_array=Sys_vals_all
            files_array=Sys_files_all
            HPs_perz = np.zeros([nthbins, 3])
        else:
            pdf_Sys_bias_perz = np.zeros([len(Sys_vals),nthbins])
            sys_array=Sys_vals
            files_array=Sys_files
            HPs_perz = np.zeros([nthbins, 1]) 
            
        indir_prefix = 'MRres%s_%s_100Sqdeg_%s_%s_%sGpAM_z%s_' %(MRres,Sys_Tag0,noise[i],Mask_Sys,Survey,Survey)
        if j != i:
            ZBlabel = 'ZBcut%s_X_ZBcut%s' %(ZBcut[i],ZBcut[j])
        else:
            ZBlabel = 'ZBcut%s' %(ZBcut[i])
        indir_Sys = '%s%s_Cosmolfid' %(indir_prefix,ZBlabel)
        filename_Sys = indir_Sys + '/%s_%s.%sGpAM.LOSAll.SS%s.SNRPDF_%sbins.bias%s-additive.dat' %(noise[i],name_Sys,
                                                                                                   Survey,SS,nthbins,Sys_Tag0)
        # ---------------------------------- Sys READING ------------------------------------------
        for s in range(len(sys_array)):
            filename_Sys2 = filename_Sys.replace(Sys_Tag0, files_array[s])
            
            if Sys == "IA":
                tmp_bias = 'IAbias%s'%Sys_vals[s]
            elif Sys == "dz" or Sys == "BaryON":
                tmp_bias = 'bias%s'%files_array[s]

            # load the pure bias:
            SNR, pdf_Sys_bias_perz[s,:] = np.loadtxt( filename_Sys2, usecols=(0,1),unpack=True)

        if Sys=="dz" and i!=j:
            # fit a 2D GP as func of dz1 & dz2:
            for t in range(nthbins):
                # rescale training points to mean 0 and stdev 1:
                pred_tmp = pdf_Sys_bias_perz[:,t].reshape(-1,1)
                scaler = preprocessing.StandardScaler().fit( pred_tmp ) # gets mean & stdev
                inpred = scaler.transform( pred_tmp )  # does transform
                #inpred = SNR[t]*pdf_Sys_bias_perz[:,t]
                GPR_Class = GPR_Emu(Sys_vals_all, inpred, np.zeros_like(inpred), Sys_vals_all)
                _,_,HPs_perz[t,:] = GPR_Class.GPRsk(np.ones(3), None, 100) # 100&500 give same answer  
        else:
            # fit a linear relation for each z&SNR bin 
            for t in range(nthbins):
                HPs_perz[t,:], cov = curve_fit(linear(x0=Sys_vals[0]), Sys_vals, SNR[t]*pdf_Sys_bias_perz[:,t])
                
        # In Make_SNR-PDF, the dz bias was computed relative to dz=-0.89; so bias=0 for this point.
        # We need to recalibrate this, so bias=0 for dz=0, as it is in the mocks & data.
        # We can't use the simple SLICS mocks for a dz=0 pred; tried this, it's outside of the dz2-5 results (for some reason).
        # So we now use the fits to compute the bias @ dz=0, and subtract this off, saving new
        # pure biases to be used in the MCMC, which will be centred on 0.
        # WHAT IS MORE... we need to read in COMBOs of dz params, since these can differ for different tomo bins
        if Sys=="dz":
            b0 = np.zeros(nthbins)
            # if it's an auto-bin, keep the simple linear fit to get the bias at dz=0
            if i==j:
                for t in range(nthbins):
                    b0[t] = linear(x0=Sys_vals[0])(0., HPs_perz[t,:]) / SNR[t]
                    # bias at dz=0 given by old fit^ (units of PDF)
                    # no need to recompute the HPs because grad is unchanged by zero-point

            # ...if X-bin, get b0 bias from 2D GPemu (dz(bin1) VS dz(bin2)) & recompute HPs
            else:
                for t in range(nthbins):
                    pred_tmp = pdf_Sys_bias_perz[:,t].reshape(-1,1)
                    # perform scaling again (mean0, stdev1):
                    scaler = preprocessing.StandardScaler().fit( pred_tmp )
                    inpred = scaler.transform( pred_tmp )      #inpred = SNR[t]*pdf_Sys_bias_perz[:,t]
                    GPR_Class = GPR_Emu(Sys_vals_all, inpred, np.zeros_like(inpred), [[0,0]])
                    b0_tmp,_,_ = GPR_Class.GPRsk(HPs_perz[t,:], None, 0) # 0 restarts (it's pre-trained)
                    b0[t] = scaler.inverse_transform( b0_tmp.reshape(-1,1) ).flatten()

                    # now recompute the HPs with the preds re-defined as centred on zero
                    innodes = Sys_vals_all #np.vstack((Sys_vals_all,[0.,0.])) # add [0,0] point to the gang
                    scaler = preprocessing.StandardScaler().fit( pred_tmp-b0[t] )
                    inpred = scaler.transform( pred_tmp-b0[t] )
                    #inpred = np.append(inpred, 0.)
                    GPR_Class = GPR_Emu(innodes, inpred, np.zeros_like(inpred), [[0,0]])
                    _,_,HPs_perz[t,:] = GPR_Class.GPRsk(HPs_perz[t,:], None, 100)
            
            for s in range(len(sys_array)):
                # subtract bias(dz=0) from previous pure biases:
                new_bias = pdf_Sys_bias_perz[s,:] - b0
                # save the new bias:                                                                                      
                new_fname_bias = filename_Sys.replace(Sys_Tag0, files_array[s]).replace('additive.dat',
                                                                                      'additive.0centred.dat')
                print(new_fname_bias)
                print("--------")
                #np.savetxt(new_fname_bias, np.c_[SNR, new_bias], header='# SNR, dz-bias[rel. to dz=0]')
                # replace the bias:
                pdf_Sys_bias_perz[s,:] = new_bias

        # store the PDF biases and HPs for this zbin:
        pdf_Sys_bias.append( pdf_Sys_bias_perz )
        HPs.append( HPs_perz ) 
        k+=1


def Ladder_Plot(x_array,
                y_array, y_array_err,
                x_array_fit, y_array_fit, 
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
                    ax1.plot( x_array_fit, y_array_fit[d,t,:], linestyle='-', color=s_m.to_rgba(SNR[t]) )
                    ax1.errorbar( x_array, y_array[d,t,:], yerr=0., #y_array_err[d,t],
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

y_limits = [-7., 6.]

# assemble the things you want to plot:
# for IA & baryons, this is just the mock points & linear model fits as a func of A_IA / b_bary
# for dz auto bins, it is this^ as well.
# But for dz X-bins, it's the mock bias & GPemu model evaluated on diag, where: dz1=dz2
if Sys=="dz":
    Sys_vals_p0 = np.insert(Sys_vals, 2, 0.) # Sys_vals but with dz=0 inserted
    GP_Sys_vals = np.column_stack((Sys_vals_p0,Sys_vals_p0)) # only used for dz GPemu in X-bins                                                                  
else:
    Sys_vals_p0 = Sys_vals

    
plot_fit = np.zeros([nzbins, nthbins, len(Sys_vals_p0)]) # add 1; show fit at 0 as well
plot_mock = np.zeros([nzbins, nthbins, len(Sys_vals)])

k=0
for i in range(len(noise)):
    for j in range(i, len(noise)):

        if Sys=="dz" and i!=j:
            # Use the 2D GP:
            for t in range(nthbins): 
                innodes = Sys_vals_all 
                # no need to add [0,0] to train set: bias(0,0)=0 behaviour recovered
                pred_tmp = pdf_Sys_bias[k][:,t].reshape(-1,1)
                scaler = preprocessing.StandardScaler().fit( pred_tmp )
                inpred = scaler.transform( pred_tmp ) 
                GPR_Class = GPR_Emu(innodes, inpred, np.zeros_like(inpred), GP_Sys_vals)
                outpred,_,_ = GPR_Class.GPRsk(HPs[k][t,:], None, 0) # 0 restarts (pre-trained)
                plot_fit[k,t,:] = scaler.inverse_transform( outpred.reshape(-1,1) ).flatten()

            
            # now pull out the mock measurements where dz1=dz2 for this X-bin (want to plot these):
            scount=0
            for s in range(len(Sys_files_all)):
                if '_X_' not in Sys_files_all[s]: # then dz1=dz2
                    plot_mock[k,:,scount] = pdf_Sys_bias[k][s,:] #SNR * pdf_Sys_bias[k][s,:]
                    #print("I: filled zbin %s of plot_mock" %k)
                    scount+=1
           
        else:
            # it's simple, just do the linear model fit:
            for t in range(nthbins):
                plot_fit[k,t,:] = linear(x0=0.)(Sys_vals_p0, HPs[k][t,:]) / SNR[t]
                plot_mock[k,t,:] = pdf_Sys_bias[k][:,t] #SNR[t]*pdf_Sys_bias[k][:,t]
                #print("II: filled zbin %s of plot_mock" %k)
                
        k+=1 # increment zbin
    

Ladder_Plot(Sys_vals,
            plot_mock, 0.,
            Sys_vals_p0, plot_fit,
            len(noise), leg_Label,
            r'$\Delta {\rm PDF}$ | $\sigma_{\rm s}=$%s' %(ss_label),
            None, y_limits, None)


