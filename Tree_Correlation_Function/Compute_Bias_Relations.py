# 03/12/2024, B. Giblin
# This reads in the sys. biases for dz or baryons, made by Plot_xipm_Ladder.py, 
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
Sys = "IA"
if Sys=="IA":
    Sys_vals = np.arange(-6., 7.0) # actual vals to use for colour scale
    Sys_labs = Sys_vals # what they'll be labelled as in legend
    Sys_files= ['IA-6.0','IA-5.0','IA-4.0','IA-3.0','IA-2.0','IA-1.0','IA0.0',
                'IA1.0','IA2.0','IA3.0','IA4.0','IA5.0','IA6.0'] #np.arange(-6., 7.0) # differences in filename
    leg_Label = r'$A_{\rm IA}$'
    c_offset = 2 #2 # helps to get colour of Sys lines visible
    sigma = np.ones(5) # dummy; means Sys_vals same for every tomo bin
     
    NLOS_Sys = 50   
    Sys_Tag0 = 'IA0.0'
    Sys_idx0 = 6 # the index in Sys_vals corresponding to SysTag0
    Mask_Sys = 'NoMask'                                                        
    name_Sys = 'test'

elif Sys=="dz":
    Sys_vals = np.array([-0.8980,-0.2293,0.2623,0.8254]) # shifts in units of sigma_z
    Sys_labs = [-0.90,-0.23,0.26,0.83]
    sigma = np.array([0.0096,0.0114,0.0116,0.0084,0.0097]) # sigma_z per tomo bin
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
    Sys_idx0 = 0 # the index in Sys_vals corresponding to SysTag0
    Mask_Sys = 'Mosaic'
    name_Sys = 'Mosaic'

elif Sys=="BaryON":
    Sys_vals = np.array([0, 1])
    Sys_labs = ["DM-only", "Baryons",]
    Sys_files= ["BaryOFF", "BaryON"]
    leg_Label = 'Baryons'
    c_offset = 0.5
    sigma = np.ones(5) # dummy; means Sys_vals same for every tomo bin
                    
    NLOS_Sys = 180
    Sys_Tag0 = 'BaryOFF'
    Sys_idx0 = 0 # the index in Sys_vals corresponding to SysTag0
    Mask_Sys = 'Mosaic'
    name_Sys = 'Mosaic'

pm = '-'
if pm == '+':
    readcol = 3
elif pm == '-':
    readcol = 4

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
    
nthbins = 9
nzbins = 15

xi_uc_Sys_bias = [] #np.zeros([ len(Sys_vals), nzbins, nthbins ]) # the pure additive bias (bias_x-bias_0)
xi_c_Sys_bias = [] #np.zeros([ len(Sys_vals), nzbins, nthbins ]) # same but for clipped


# define function for linear fit (just grad; defined so bias=0 at x[0])
from scipy.optimize import curve_fit

# linear model used to fit sys-bias per bin (e.g. per theta for 2pcf) as a func of Sys_Nodes x (e.g. dz).
# Because curve_fit doesnt play nice with fixed params (x0), need an outer & an inner loop to provide this.
def linear(x0):
    def inner_linear(x, m):
        return m*(x-x0)
    return inner_linear

def Apply_Norm(pred):
    # rescale training points to mean 0 and stdev 1:
    scaler = preprocessing.StandardScaler().fit( pred )
    scaled_pred = scaler.transform( pred )  
    return scaler, scaled_pred

# store the gradient (hyperparam) & its error for each z&theta bin
HPs_uc = []     # unclipped
HPs_c = []      # clipped
sys_array_store = []
k=0
for i in range(len(noise)):
    for j in range(i, len(noise)):
        if Sys=="dz" and i!=j:
            xi_uc_Sys_bias_perz = np.zeros([len(Sys_vals_all),nthbins])
            xi_c_Sys_bias_perz = np.zeros([len(Sys_vals_all),nthbins])
            HPs_uc_perz = np.zeros([nthbins, 3])
            HPs_c_perz = np.zeros([nthbins, 3])
            sys_array=Sys_vals_all * [sigma[i],sigma[j]] # convert from /sigma units to raw 
            files_array=Sys_files_all
        else:
            xi_c_Sys_bias_perz = np.zeros([len(Sys_vals),nthbins])
            xi_uc_Sys_bias_perz = np.zeros([len(Sys_vals),nthbins])
            HPs_c_perz = np.zeros([nthbins, 1])
            HPs_uc_perz = np.zeros([nthbins, 1])
            sys_array=Sys_vals * sigma[i]
            files_array=Sys_files
        sys_array_store.append(sys_array)

        indir_prefix = 'MRres%s_%s_100Sqdeg_%s_%s_%sGpAM_z%s_' %(MRres,Sys_Tag0,noise[i],Mask_Sys,Survey,Survey)
        if j != i:
            ZBlabel = 'ZBcut%s_X_ZBcut%s' %(ZBcut[i],ZBcut[j])
        else:
            ZBlabel = 'ZBcut%s' %(ZBcut[i])
        indir_Sys = '%s%s_Cosmolfid/ThBins%s/NLOS%s' %(indir_prefix,ZBlabel,nthbins,NLOS_Sys)
        # unclipped & clipped:
        filename_Sys = indir_Sys + '/%s_%s.%sGpAM.NLOS%s.AverageCF%s.asc' %(noise[i],name_Sys,Survey,NLOS_Sys,pm)
        filename_c_Sys = filename_Sys.replace( 'AverageCF%s' %pm, 'SS%s.rCLIP_X3sigma.AverageCF%s'%(SS,pm) )

        # ---------------------------------- Sys READING ------------------------------------------
        for s in range(len(sys_array)):
            # unclipped:
            filename_Sys2 = filename_Sys.replace(Sys_Tag0, files_array[s])
            # clipped:
            filename_c_Sys2 = filename_c_Sys.replace(Sys_Tag0, files_array[s])

            # now load the pure bias:
            tmp_bias = 'bias%s'%files_array[s]
            fname_bias = filename_Sys2.replace('.asc','.%s-additive.asc'%tmp_bias)
            fname_bias_c = filename_c_Sys2.replace('.asc','.%s-additive.asc'%tmp_bias)
            theta, xi_uc_Sys_bias_perz[s,:] = np.loadtxt( fname_bias, usecols=(0,1),unpack=True)
            xi_c_Sys_bias_perz[s,:] = np.loadtxt( fname_bias_c, usecols=(1,),unpack=True)

        for t in range(nthbins):
            inpred_uc = 1e4*theta[t]*xi_uc_Sys_bias_perz[:,t].reshape(-1,1)
            inpred_c = 1e4*theta[t]*xi_c_Sys_bias_perz[:,t].reshape(-1,1)

            if Sys=="dz" and i!=j:
                # fit a 2D GP as func of dz1 & dz2:
                GPR_Class = GPR_Emu(sys_array, inpred_uc, np.zeros_like(inpred_uc), sys_array)
                _,_,HPs_uc_perz[t,:] = GPR_Class.GPRsk(np.ones(3), None, 100)
                GPR_Class = GPR_Emu(sys_array, inpred_c, np.zeros_like(inpred_c), sys_array)
                _,_,HPs_c_perz[t,:] = GPR_Class.GPRsk(np.ones(3), None, 100) # 100&500 give same answer

            else:
                # linear fit
                HPs_uc_perz[t,:], cov = curve_fit(linear(x0=sys_array[Sys_idx0]), sys_array, inpred_uc.flatten()) 
                HPs_c_perz[t,:], cov = curve_fit(linear(x0=sys_array[Sys_idx0]), sys_array, inpred_c.flatten()) 

        # In Plot_xipm_Ladder.py, the dz bias was computed relative to dz=-0.89; so bias=0 for this point.
        # We need to recalibrate this, so bias=0 for dz=0, as it is in the mocks & data.
        # We can't use the simple SLICS mocks for a dz=0 pred; tried this, it's outside of the dz2-5 results (for some reason).
        # So we now use the fits to compute the bias @ dz=0, and subtract this off, saving new
        # pure biases to be used in the MCMC, which will be centred on 0.
        if Sys=="dz":
            b0 = np.zeros(nthbins)
            b0_c = np.zeros(nthbins)

            for t in range(nthbins):
                inpred_uc = 1e4*theta[t]*xi_uc_Sys_bias_perz[:,t].reshape(-1,1)
                inpred_c = 1e4*theta[t]*xi_c_Sys_bias_perz[:,t].reshape(-1,1)

                if i==j:
                    # an auto-bin, keep the simple linear fit to get the bias at dz=0
                    b0[t] = linear(x0=sys_array[Sys_idx0])(0., HPs_uc_perz[t,:]) / (1e4*theta[t])
                    b0_c[t] = linear(x0=sys_array[Sys_idx0])(0., HPs_c_perz[t,:]) / (1e4*theta[t])
                    # shift training predictions so centred on 0:
                    inpred_uc = 1e4*theta[t]*(xi_uc_Sys_bias_perz[:,t] - b0[t]).reshape(-1,1)
                    inpred_c = 1e4*theta[t]*(xi_c_Sys_bias_perz[:,t] - b0_c[t]).reshape(-1,1)
                    # retrain (may not be change HPs, but being careful):
                    HPs_uc_perz[t,:], cov = curve_fit(linear(x0=0.), sys_array, inpred_uc.flatten())
                    HPs_c_perz[t,:], cov = curve_fit(linear(x0=0.), sys_array, inpred_c.flatten()) # HPs now for b0-corrected PDFs
                    
                else:
                    GPR_Class = GPR_Emu(sys_array, inpred_uc, np.zeros_like(inpred_uc), [[0,0]])
                    b0[t] = GPR_Class.GPRsk(HPs_uc_perz[t,:], None, 0)[0][0] / (1e4*theta[t]) # 0 restarts (it's pre-trained)
                    GPR_Class = GPR_Emu(sys_array, inpred_c, np.zeros_like(inpred_c), [[0,0]])
                    b0_c[t] = GPR_Class.GPRsk(HPs_c_perz[t,:], None, 0)[0][0] / (1e4*theta[t])
                    # shift training predictions so centred on 0:
                    inpred_uc = 1e4*theta[t]*(xi_uc_Sys_bias_perz[:,t] - b0[t]).reshape(-1,1)
                    inpred_c = 1e4*theta[t]*(xi_c_Sys_bias_perz[:,t] - b0_c[t]).reshape(-1,1)
                    # retrain:
                    GPR_Class = GPR_Emu(sys_array, inpred_uc, np.zeros_like(inpred_uc), [[0,0]])
                    _,_,HPs_uc_perz[t,:] = GPR_Class.GPRsk(HPs_uc_perz[t,:], None, 100)
                    GPR_Class = GPR_Emu(sys_array, inpred_c, np.zeros_like(inpred_c), [[0,0]])
                    _,_,HPs_c_perz[t,:] = GPR_Class.GPRsk(HPs_c_perz[t,:], None, 100)


            for s in range(len(sys_array)):
                # subtract bias(dz=0) from previous pure biases:
                new_bias = xi_uc_Sys_bias_perz[s,:] - b0
                new_bias_c = xi_c_Sys_bias_perz[s,:] - b0_c
                # save the new bias:
                new_fname_bias = fname_bias.replace(files_array[-1], files_array[s]).replace('additive.asc',
                                                                                          'additive.0centred.asc')
                new_fname_bias_c = fname_bias_c.replace(files_array[-1], files_array[s]).replace('additive.asc',
                                                                                              'additive.0centred.asc')
                np.savetxt(new_fname_bias, np.c_[theta, new_bias], header='# theta[arcmin], dz-bias[rel. to dz=0]')
                np.savetxt(new_fname_bias_c, np.c_[theta, new_bias_c], header='# theta[arcmin], dz-bias[rel. to dz=0]')
                
                print(new_fname_bias)
                #print(new_fname_bias_c)
                #print("-------------------")
                # replace biases:
                xi_uc_Sys_bias_perz[s,:] = new_bias
                xi_c_Sys_bias_perz[s,:] = new_bias_c

        xi_uc_Sys_bias.append( xi_uc_Sys_bias_perz )
        xi_c_Sys_bias.append( xi_c_Sys_bias_perz )
        HPs_uc.append(HPs_uc_perz)
        HPs_c.append(HPs_c_perz)
        k+=1

# Load the SLICS err corresponding to a complete (shuffled) survey footprint
# this will be used as the error on the Sys measurement:
indir_S = 'MRres%s_100Sqdeg_SNCycle_%s_%sGpAM_z%s_ZBcut0.1-0.3_X_ZBcut0.1-0.3/ThBins%s/NLOS%s/' %(MRres,Mask,Survey,Survey,
                                                                                                  nthbins,NLOS)
# unclipped:
filename_cov = '%s_%s.%sGpAM.NLOS%s-Shuffle0.ORIG.UCxUC%s.CovMat.zbins11_12_13_14_15_22_23_24_25_33_34_35_44_45_55.npy' %(noise[i],filetag,Survey,NLOS,pm)
filename_cov = indir_S + filename_cov
xi_err = np.sqrt( np.diag( np.load(filename_cov) ) ) # stdev for one K1000 survey                               
xi_err = np.reshape(xi_err, (nzbins,nthbins))        # reshape
xi_err = xi_err / np.sqrt(NLOS_Sys/18.)              # make it error on the mean of 1 survey
# clipped:
filename_c_cov = filename_cov.replace( 'UCxUC%s', 'SS%s.rCLIP_X3sigma.CxC%s'%(SS,pm) )
xi_c_err = np.sqrt( np.diag( np.load(filename_c_cov) ) )
xi_c_err = np.reshape(xi_c_err, (nzbins,nthbins))
xi_c_err = xi_c_err / np.sqrt(NLOS_Sys/18.)


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
                for t in [0,3,6,8]: #range(nthbins):
                    #ax1.plot( x_array, y_array[:,d,t], color=s_m.to_rgba(np.log10(theta[t])) )
                    ax1.plot( x_array_fit, y_array_fit[d,t,:], linestyle='-', color=s_m.to_rgba(np.log10(theta[t]) ) )
                    ax1.errorbar( x_array, y_array[d,t,:], yerr=0., #y_array_err[d,t],
                                  fmt='o', ecolor=s_m.to_rgba(np.log10(theta[t]) ),
                                  markerfacecolor=s_m.to_rgba(np.log10(theta[t]) ),
                                  markeredgecolor=s_m.to_rgba(np.log10(theta[t]) ) )
                    
                ax1.text( 0.75, 0.8, r'%s-%s' %(j+1,i+1),
                                 horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes)
                #ax1.set_yscale('log')
                #ax1.set_xscale('log')
                linthresh=0.1
                #ax1.set_yscale('symlog', linthresh=linthresh, linscale=1.0)
                ax1.plot([(x_array*sigma[i]).min(),(x_array*sigma[i]).max()], [-1*linthresh,-1*linthresh], 'k:')
                ax1.plot([(x_array*sigma[i]).min(),(x_array*sigma[i]).max()], [linthresh,linthresh], 'k:')
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
                elif Sys != "dz":
                    # dont kill xvals if working with dz, since x vals differ per tomo bin.
                    ax1.set_xticks([])
                d+=1

    if Sys != "dz":
        fig.subplots_adjust(hspace=0, wspace=0)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(s_m, cax=cbar_ax, label=r'$\log_{10} \theta \, [\rm{arcmin}]$')
    plt.savefig(savename)
    plt.show()
    return

# also make a colour scale for theta
c_m = matplotlib.cm.jet
norm = matplotlib.colors.Normalize(vmin=np.log10(theta).min(), vmax=np.log10(theta).max()) 
s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
s_m.set_array([])


if pm == '+':
    if Sys=="dz":
        y_limits=[-0.09, 0.12]
    elif Sys=="IA":
        y_limits=[-0.49, 0.5]
    elif Sys=="BaryON":
        y_limits=[-0.04, 0.009]

elif pm == '-':
    if Sys=="dz":
        y_limits=[-0.04, 0.06]
    elif Sys=="IA":
        y_limits=[-0.24, 0.26]
    elif Sys=="BaryON":
        y_limits=[-0.04, 0.009]
    

# assemble the things you want to plot:
# for IA & baryons, this is just the mock points & linear model fits as a func of A_IA / b_bary
# for dz auto bins, it is this^ as well.
# But for dz X-bins, it's the mock bias & GPemu model evaluated on diag, where: dz1=dz2
if Sys=="dz":
    Sys_vals_p0 = np.insert(Sys_vals, 2, 0.) # Sys_vals but with dz=0 inserted
    GP_Sys_vals = np.column_stack((Sys_vals_p0,Sys_vals_p0)) # only used for dz GPemu in X-bins                                                                  
else:
    Sys_vals_p0 = Sys_vals

plot_fit_uc = np.zeros([nzbins, nthbins, len(Sys_vals_p0)]) # adding 1 element in dz case, to show fit at 0 as well
plot_mock_uc = np.zeros([nzbins, nthbins, len(Sys_vals)])
plot_fit_c = np.zeros_like( plot_fit_uc )
plot_mock_c = np.zeros_like( plot_mock_uc )

k=0
for i in range(len(noise)):
    for j in range(i, len(noise)):
        
        if Sys=="dz" and i!=j:
            # Use the 2D GP:
            for t in range(nthbins): 
                # no need to add [0,0] to train set: bias(0,0)=0 behaviour recovered
                inpred_uc = 1e4*theta[t]*xi_uc_Sys_bias[k][:,t].reshape(-1,1)
                inpred_c = 1e4*theta[t]*xi_c_Sys_bias[k][:,t].reshape(-1,1)
                # unclipped:
                GPR_Class = GPR_Emu(sys_array_store[k], inpred_uc, np.zeros_like(inpred_uc), GP_Sys_vals*[sigma[i],sigma[j]] )
                plot_fit_uc[k,t,:] = GPR_Class.GPRsk(HPs_uc[k][t,:], None, 0)[0].flatten() 
                # clipped:
                GPR_Class = GPR_Emu(sys_array_store[k], inpred_c, np.zeros_like(inpred_c), GP_Sys_vals*[sigma[i],sigma[j]] )
                plot_fit_c[k,t,:] = GPR_Class.GPRsk(HPs_c[k][t,:], None, 0)[0].flatten() 
            
            # now pull out the mock measurements where dz1=dz2 for this X-bin (want to plot these):
            scount=0
            for s in range(len(Sys_files_all)):
                if '_X_' not in Sys_files_all[s]: # then dz1=dz2
                    plot_mock_uc[k,:,scount] = xi_uc_Sys_bias[k][s,:] * 1e4 * theta
                    plot_mock_c[k,:,scount] = xi_c_Sys_bias[k][s,:] * 1e4 * theta
                    #print("I: filled zbin %s of plot_mock" %k)
                    scount+=1
           
        else:
            # it's simple, just do the linear model fit:
            for t in range(nthbins):
                plot_fit_uc[k,t,:] = linear(x0=0.)(Sys_vals_p0*sigma[i], HPs_uc[k][t,:]) 
                plot_fit_c[k,t,:] = linear(x0=0.)(Sys_vals_p0*sigma[i], HPs_c[k][t,:]) 
                plot_mock_uc[k,t,:] = xi_uc_Sys_bias[k][:,t] * 1e4 * theta[t]
                plot_mock_c[k,t,:] = xi_c_Sys_bias[k][:,t] * 1e4 * theta[t]
                #print("II: filled zbin %s of plot_mock" %k)
                
        k+=1 # increment zbin

# unclipped
Ladder_Plot(Sys_vals,
            plot_mock_uc, 0.,
            Sys_vals_p0, plot_fit_uc,
            len(noise), leg_Label,
            r'$\theta \times \Delta \xi_{%s}^{\rm uc} \, [10^{-4} \rm{arcmin}]$' %pm,
            None, y_limits, '/home/bengib/unclipped.png')

# clipped
Ladder_Plot(Sys_vals,
            plot_mock_c, 0.,
            Sys_vals_p0, plot_fit_c,
            len(noise), leg_Label,
            r'$\theta \times \Delta \xi_{%s}^{\rm c} \, [10^{-4} \rm{arcmin}]$' %pm,
            None, y_limits, '/home/bengib/clipped.png')
