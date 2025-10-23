# 03/08/2021, B. Giblin
# Plot the PDF(SNR) for all redshift bin combinations as a ladder plot

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

Plot_Ratio = True      # plot the ratio of the measurements to the fid. cosmology.
Plot_Sys = False          # if True, sets cbar to Sys values (IA/dz/Bary; overwriting S8/Delta_S8 cbar)
Sys = "IA"
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
    Sys_vals = [-0.8980,-0.2293,0.2623,0.8254]
    Sys_labs = [-0.90,-0.23,0.26,0.83]
    Sys_files= ['dz5','dz3','dz4','dz2']
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
    Sys_vals = [0, 1]
    Sys_labs = ["DM-only", "Baryons",]
    Sys_files= ["BaryOFF", "BaryON"]
    leg_Label = 'Baryons'
    c_offset = 0.5
                    
    NLOS_Sys = 180
    Sys_Tag0 = 'BaryOFF'
    Mask_Sys = 'Mosaic'
    name_Sys = 'Mosaic'

pm = '-'
if pm == '+':
    readcol = 3
elif pm == '-':
    readcol = 4

cosmol = []
for i in range(25):
    cosmol.append('%s' %i)
cosmol.append('fid')
#cosmol = 'fid'

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
    NLOS = 3906 #715
elif Survey == 'LSST':
    noise = [ 'SN0.28', 'SN0.28', 'SN0.28', 'SN0.28', 'SN0.28']
    NLOS = 616
ZBcut = ['0.1-0.3', '0.3-0.5', '0.5-0.7', '0.7-0.9', '0.9-1.2'] #['0.1-1.2']
    
nthbins = 9
nzbins = 15
xi_uc = np.zeros([ len(cosmol), nzbins, nthbins ])         # cosmoSLICS mean
theta_CS = np.zeros([ len(cosmol), nzbins, nthbins ])

xi_c = np.zeros_like( xi_uc )                           # clipped cosmoSLICS mean

xi_uc_S  = np.zeros([ nzbins, nthbins ])                      # SLICS mean
theta_S = np.zeros([ nzbins, nthbins ])
xi_c_S = np.zeros_like( xi_uc_S )                             # clipped SLICS mean

xi_uc_all = np.zeros([ len(cosmol), 50, nzbins, nthbins ]) # cosmoSLICS all realisations
theta_CS_all = np.zeros([ len(cosmol), 50, nzbins, nthbins ])

xi_uc_S_all = np.zeros([ NLOS, nzbins, nthbins ])              # SLICS all realisations
theta_S_all = np.zeros([ NLOS, nzbins, nthbins ])

# Sys (IA/dz/Baryon) measurements!
xi_uc_Sys = np.zeros([ nzbins, nthbins ])           # IA0.0/dz0/BaryOFF measurement 
xi_c_Sys = np.zeros([ nzbins, nthbins ])            # same but clipped:
xi_uc_Sys2 = []         # All other sys measurements
xi_c_Sys2 = [] 
xi_uc_Sys_bias = []     # the pure additive bias (bias_x-bias_0)
xi_c_Sys_bias = []      # same but for clipped
# CS-fid CONTAMINATED by Sys
xi_uc_Sys_add = []      # CS-fid with Sys manually added in (additive)
xi_c_Sys_add = []       # same, but for clipped:
# SLICS CONTAMINATED by Sys
xi_uc_S_Sys_add = [] 
xi_c_S_Sys_add = [] 

xi_err  = np.zeros([ nzbins, nthbins ])                    # SLICS error
xi_c_err = np.zeros_like( xi_err )                         # clipped SLICS error

# KiDS
xi_uc_K = np.zeros_like( xi_uc_S )
xi_c_K = np.zeros_like( xi_uc_S )

for c in range(len(cosmol)):
    k=0
    print("On cosmology ", cosmol[c])
    for i in range(len(noise)):
        for j in range(i, len(noise)):
            indir_prefix = 'MRres%s_100Sqdeg_%s_%s_%sGpAM_z%s_ZBcut%s' %(MRres,noise[i],Mask,Survey,Survey,ZBcut[i])
            if j != i:
                indir_prefix += '_X_ZBcut%s' %ZBcut[j]
            # Load cosmoSLICS
            indir_CS = indir_prefix + '_Cosmol%s/ThBins%s' %(cosmol[c], nthbins)
            
            # find all realisations
            all_filenames = glob.glob('%s/%s_%s.%sGpAM.LOS*.ORIG.CorrFun.asc' %(indir_CS,noise[i],filetag,Survey))
            #for l in range(len(all_filenames)):
            #    theta_CS_all[c,l,k,:], xi_uc_all[c,l,k,:] = np.loadtxt( all_filenames[l], usecols=(1,readcol), unpack=True )
            
            indir_CS += '/NLOS%s' %NLOS_CS     # cosmoSLICS
            filename_CS = indir_CS + '/%s_%s.%sGpAM.NLOS%s.AverageCF%s.asc' %(noise[i],filetag,Survey,NLOS_CS,pm)
            #print( filename_CS )
            theta_CS[c,k,:], xi_uc[c,k,:] = np.loadtxt(filename_CS, usecols=(0,1), unpack=True)
            filename_c_CS = filename_CS.replace( 'AverageCF%s' %pm, 'SS%s.rCLIP_X3sigma.AverageCF%s'%(SS,pm) )
            xi_c[c,k,:] = np.loadtxt( filename_c_CS, usecols=(1,), unpack=True )
            
            if c == len(cosmol)-1:
                print("-----------------------------------  ZBcut%s_X_%s --------------------------------- " %(ZBcut[i],ZBcut[j]))
                # Load KiDS:
                filename_KiDS = filename_CS.replace('fid','KiDS1000').replace('NLOS900','NLOS18').replace(noise[i],'SNKiDS')
                filename_c_KiDS = filename_c_CS.replace('fid','KiDS1000').replace('NLOS900','NLOS18').replace(noise[i],'SNKiDS')
                xi_uc_K[k,:]   = np.loadtxt( filename_KiDS, usecols=(1,), unpack=True )
                xi_c_K[k,:] = np.loadtxt( filename_c_KiDS, usecols=(1,), unpack=True )
                
                # Last time load SLICS
                indir_S  = indir_prefix + '/ThBins%s' %nthbins                          # SLICS

                # find all realisations
                all_filenames = glob.glob('%s/%s_%s.%sGpAM.LOS*.ORIG.CorrFun.asc' %(indir_S,noise[i],filetag,Survey))
                #for l in range(len(all_filenames)):
                #    theta_S_all[l,k,:], xi_uc_S_all[l,k,:] = np.loadtxt( all_filenames[l], usecols=(1,readcol), unpack=True )

                indir_S += '/NLOS%s' %NLOS
                filename_S  = indir_S + '/%s_%s.%sGpAM.NLOS%s.AverageCF%s.asc' %(noise[i],filetag,Survey,NLOS,pm)
                theta_S[k,:], xi_uc_S[k,:] = np.loadtxt(filename_S, usecols=(0,1), unpack=True)
                # same for clipped:
                filename_c_S = filename_S.replace( 'AverageCF%s' %pm, 'SS%s.rCLIP_X3sigma.AverageCF%s'%(SS,pm) )
                xi_c_S[k,:] = np.loadtxt( filename_c_S, usecols=(1,), unpack=True )
                
                filename_cov = indir_S + '/%s_%s.%sGpAM.NLOS%s.UCxUC%s.CovMat.npy' %(noise[i],filetag,Survey,NLOS,pm)
                xi_err[k,:] = np.sqrt( np.diag( np.load(filename_cov) ) / 10. ) # scaled from 100 --> 1000 deg^2
                filename_c_cov = filename_cov.replace( 'UCxUC%s', 'SS%s.rCLIP_X3sigma.CxC%s'%(SS,pm) )
                xi_c_err[k,:] = np.sqrt( np.diag( np.load(filename_c_cov) ) / 10. )
                
                # Check you're reading in the correct xi_pm!
                #print( " --- cosmoSLICS --- " )
                #print( indir_CS )

                #print( " --- SLICS ---" )
                #print( indir_S )

                # ---------------------------------- Sys READING ------------------------------------------

                # dz X-bins have MORE bias measurements than auto-bins! SO we need to change array size:
                if Sys=="dz" and i!=j:
                    files_array=Sys_files_all
                else:
                    files_array=Sys_files

                # actual sys measurements:
                xi_uc_Sys2_perz = np.zeros([len(files_array),nthbins]) 
                xi_c_Sys2_perz = np.zeros([len(files_array),nthbins]) 
                # the pure additive bias (bias_x-bias_0)
                xi_uc_Sys_bias_perz = np.zeros([len(files_array),nthbins])  
                xi_c_Sys_bias_perz = np.zeros([len(files_array),nthbins]) # clipped
                # cosmoSLICS CONTAMINATED by Sys
                xi_uc_Sys_add_perz = np.zeros([len(files_array),nthbins]) 
                xi_c_Sys_add_perz = np.zeros([len(files_array),nthbins]) 
                # SLICS CONTAMINATED by Sys
                xi_uc_S_Sys_add_perz = np.zeros([len(files_array),nthbins])
                xi_c_S_Sys_add_perz = np.zeros([len(files_array),nthbins])

                indir_prefix = 'MRres%s_%s_100Sqdeg_%s_%s_%sGpAM_z%s_' %(MRres,Sys_Tag0,noise[i],Mask_Sys,Survey,Survey)
                if j != i:
                    ZBlabel = 'ZBcut%s_X_ZBcut%s' %(ZBcut[i],ZBcut[j])
                else:
                    ZBlabel = 'ZBcut%s' %(ZBcut[i])
                indir_Sys = '%s%s_Cosmolfid/ThBins%s/NLOS%s' %(indir_prefix,ZBlabel,nthbins,NLOS_Sys)
                # unclipped & clipped:
                filename_Sys = indir_Sys + '/%s_%s.%sGpAM.NLOS%s.AverageCF%s.asc' %(noise[i],name_Sys,Survey,NLOS_Sys,pm)
                filename_c_Sys = filename_Sys.replace( 'AverageCF%s' %pm, 'SS%s.rCLIP_X3sigma.AverageCF%s'%(SS,pm) )

                # load benchmark Sys (IA0.0 / dz0 / BaryOFF) unclipped & clipped
                xi_uc_Sys[k,:] = np.loadtxt(filename_Sys, usecols=(1,), unpack=True)
                xi_c_Sys[k,:] = np.loadtxt(filename_c_Sys, usecols=(1,), unpack=True) 
                
                for s in range(len(files_array)):
                    # unclipped:
                    filename_Sys2 = filename_Sys.replace(Sys_Tag0, files_array[s])
                    xi_uc_Sys2_perz[s,:] = np.loadtxt(filename_Sys2, usecols=(1,), unpack=True)
                    # clipped:
                    filename_c_Sys2 = filename_c_Sys.replace(Sys_Tag0, files_array[s])
                    xi_c_Sys2_perz[s,:] = np.loadtxt(filename_c_Sys2, usecols=(1,), unpack=True)

                    # compute the pure Sys bias (additive & multiplicative):
                    bias_add = (xi_uc_Sys2_perz[s,:] - xi_uc_Sys[k,:]) # unclipped 
                    bias_c_add = (xi_c_Sys2_perz[s,:] - xi_c_Sys[k,:]) #clipped
                    # store the additive bias
                    xi_uc_Sys_bias_perz[s,:] = bias_add
                    xi_c_Sys_bias_perz[s,:] = bias_c_add

                    # save the pure bias:
                    tmp_bias = 'bias%s'%files_array[s]
                    np.savetxt( filename_Sys2.replace('.asc','.%s-additive.asc'%tmp_bias), np.c_[theta_CS[c,k,:], bias_add] )
                    np.savetxt( filename_c_Sys2.replace('.asc','.%s-additive.asc'%tmp_bias),
                                np.c_[theta_CS[c,k,:], bias_c_add] )
                    print("Saving this bias files (& its clipped analog):")
                    print(filename_Sys2.replace('.asc','.%s-additive.asc'%tmp_bias))
                    #print(filename_c_Sys2.replace('.asc','.%s-additive.asc'%tmp_bias))
                    
                    # manually add/scale in the Sys contribution to CS-fid
                    # unclipped:
                    xi_uc_Sys_add_perz[s,:] = xi_uc[-1,k,:] + bias_add
                    #clipped:
                    xi_c_Sys_add_perz[s,:] = xi_c[-1,k,:] + bias_c_add
                    # same for SLICS:
                    xi_uc_S_Sys_add_perz[s,:] = xi_uc_S[k,:] + bias_add
                    xi_c_S_Sys_add_perz[s,:] = xi_c_S[k,:] + bias_c_add
                    
                    # now save them... (modify the cosmoSLICS-fid savename)
                    filename_add = filename_CS.replace('AverageCF%s.asc' %(pm), 'AverageCF%s.%s-added.asc' %(pm,tmp_bias))
                    filename_c_add = filename_c_CS.replace('AverageCF%s.asc' %(pm), 'AverageCF%s.%s-added.asc' %(pm,tmp_bias))
                    # unclipped:
                    np.savetxt( filename_add, np.c_[theta_CS[c,k,:], xi_uc_Sys_add_perz[s,:]] )
                    # clipped
                    np.savetxt( filename_c_add, np.c_[theta_CS[c,k,:], xi_c_Sys_add_perz[s,:]] )
                    # And the same for SLICS:
                    filename_add = filename_S.replace('AverageCF%s.asc' %(pm), 'AverageCF%s.%s-added.asc' %(pm,tmp_bias))
                    filename_c_add = filename_c_S.replace('AverageCF%s.asc' %(pm), 'AverageCF%s.%s-added.asc' %(pm,tmp_bias))
                    np.savetxt( filename_add, np.c_[theta_CS[c,k,:], xi_uc_S_Sys_add_perz[s,:]] )
                    np.savetxt( filename_c_add, np.c_[theta_CS[c,k,:], xi_c_S_Sys_add_perz[s,:]] )

                # store all the per-z sys arrays:
                xi_uc_Sys2.append( xi_uc_Sys2_perz )
                xi_c_Sys2.append( xi_c_Sys2_perz )
                # the pure additive bias (bias_x-bias_0)
                xi_uc_Sys_bias.append( xi_uc_Sys_bias_perz )
                xi_c_Sys_bias.append( xi_c_Sys_bias_perz ) 
                # cosmoSLICS CONTAMINATED by Sys
                xi_uc_Sys_add.append( xi_uc_Sys_add_perz )
                xi_c_Sys_add.append( xi_c_Sys_add_perz )
                # SLICS CONTAMINATED by Sys
                xi_uc_S_Sys_add.append( xi_uc_S_Sys_add_perz )
                xi_c_S_Sys_add.append( xi_c_S_Sys_add_perz ) 

            k+=1




# avg all the theta arrays
avg_theta_CS_all = np.mean( theta_CS_all, axis=1 ) # avg all realisations, per cos & tomo
                                                   # encouraging we basically get the same theta for all cos and tomo's.
avg_theta_CS = np.mean( avg_theta_CS_all, axis=(0,1) ) # overall avg for cosmoSLICS
                                                   
avg_theta_S_all = np.mean( theta_S_all, axis=0   ) # avg all realisations, per tomo
avg_theta_S = np.mean( theta_S, axis=0 )           # overall avg for SLICS - basically the same as the theta per tomo.

#theta_CS_all[c,l,k,:]
#theta_S_all[l,k,:]
#theta_S[k,:]
#theta_CS_fid_all_JHD[l,k,:]




def Ladder_Plot(x_array,
                y_array, #y_array_all,
                x_array_data,
                y_array_data, y_array_err,
                y_array_Sys,
                #y_array_data_all,
                #y_array_data,
                #y_array_data2, #y_array_data2_all,
                #y_array_data3,
                n_row_col, x_label, y_label, x_lims, y_lims, savename):

    fig = plt.figure(figsize = (13.5,10)) #figsize = (20,14)
    gs1 = gridspec.GridSpec(n_row_col, n_row_col)
    d=0     # data number
    
    axes = []
    for j in range(n_row_col):   # scroll across columns 
        for i in range(n_row_col):       # scroll down rows                                                                   

            #axes.append( plt.subplot(gs1[p]) )
            ax1 = plt.subplot(gs1[i,j])
            if j>i:
                # Dont plot panels above the diagonal                                                                        
                ax1.axis('off')
            else:
                
                # plot mean for all cosmologies
                for c in range(26): #y_array.shape[0]):
                    if Plot_Ratio:
                        #if S8[c] in [0.6101, 0.6615, 0.7232, 0.7821, 0.8321, 0.8947]:
                        ax1.plot( x_array, y_array[c,d,:], color=s_m2.to_rgba(diff_S8[c]) )
                    else:
                        ax1.plot( x_array, y_array[c,d,:], color=s_m.to_rgba(S8[c]) ) 
                if Plot_Ratio:
                    ax1.plot( x_array, np.ones_like(x_array), 'k:' )
                
                if Plot_Sys:
                    if Sys=="dz" and i!=j:
                        plot_sys_idx = [0,5,10,15] # only plot where dz_i==dz_j
                    else:
                        plot_sys_idx = range(len(Sys_vals))
                    for s in range(len(Sys_vals)):
                        idx = plot_sys_idx[s]
                        ax1.plot( x_array, y_array_Sys[d][idx,:], color=s_m3.to_rgba(Sys_vals[s]) )
                    
                # plot all realisations for final cosmology
                #for l in range(y_array_all.shape[1]):
                #    ax1.plot( x_array, y_array_all[-1,l,d,:], color='limegreen' )
                # plot mean for final cosmology
                #ax1.plot( x_array, y_array[-1,d,:], color=s_m.to_rgba(S8[-1]) )
                    
                # plot all realisations of data
                #for l in range(y_array_data2_all.shape[0]):
                #    ax1.plot( x_array, y_array_data2_all[l,d,:], color='dimgrey' )

                # plot data
                ax1.plot( x_array_data, y_array_data[d,:], color='magenta', linewidth=5 )
                #ax1.plot( x_array_data, y_array_data2[d,:], color='purple', linewidth=2 )
                #ax1.plot( x_array_data, y_array_data3[d,:], color='black', linewidth=2 )

                #ax1.errorbar( x_array_data, y_array_data3[d,:], yerr=y_array_err[d,:], color='purple', linewidth=2 )
                #ax1.errorbar( x_array_data2, y_array_data2[d,:], color='cyan', linewidth=2 )

                ax1.text( 0.75, 0.8, r'%s-%s' %(j+1,i+1),
                                 horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes)
                #ax1.set_yscale('log')
                ax1.set_xscale('log')
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
    if Plot_Sys:
        fig.colorbar(s_m3, cax=cbar_ax, label=leg_Label)
    elif Plot_Ratio:
        fig.colorbar(s_m2, cax=cbar_ax, label=r'$\Delta S_8$')
    else:
        fig.colorbar(s_m, cax=cbar_ax, label=r'$S_8$')
        
    #plt.savefig(savename)
    plt.show()
    return

# Read in the S8 values for the cosmoSLICS
S8 = np.loadtxt('/home/bengib/cosmoSLICS/Cosmologies/SLICS_Cosmol_Table_ExtraCol_25plusFid.txt',
                usecols=(1,), unpack=True)
norm = matplotlib.colors.Normalize(vmin=S8.min(), vmax=S8.max())
c_m = matplotlib.cm.jet #cool
s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
s_m.set_array([])

# also make one of these for difference in S_8 relative to fid cosmology:
diff_S8 = S8-S8[-1]
norm2 = matplotlib.colors.Normalize(vmin=diff_S8.min(), vmax=diff_S8.max())
s_m2 = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm2)
s_m2.set_array([])

# also make a colour scale for Sys
c_m3 = matplotlib.cm.Greys
norm3 = matplotlib.colors.Normalize(vmin=np.min(Sys_vals)-c_offset, vmax=np.max(Sys_vals)) #min-2 so no line's whte
s_m3 = matplotlib.cm.ScalarMappable(cmap=c_m3, norm=norm3)
s_m3.set_array([])


if pm == '+':
    y_limits_uc=[-0.9, 6.9]  #[-0.5,4.5] #[-6., 14.]
    y_limits_c=[-7.9, 7.9]
    ylog_limits=[1e-8, 1e-3] 

elif pm == '-':
    y_limits_uc=[-1.9,3.] #[-9., 9.]
    y_limits_c=[-1.9,3.]
    ylog_limits=[1e-9, 1e-4]  
    
xi_uc_fid = np.reshape(xi_uc[-1,:,:], (1, xi_uc.shape[1], xi_uc.shape[2]))



if Plot_Ratio:
    # Edit this section to vary whether you plot ratio
    # of all unclipped measurements to unclipped cosmoSLICS_fid
    # (they're calculated in the for loop below...)
    # OR you plot the ratio of clipped to uclipped,
    # which is calculated with these 3 lines here:
    Diff_c = 100*(xi_c / xi_uc -1) # ALL Diff_ == PERCENTAGE DIFFERENCES
    Diff_c_S = 100*(xi_c_S / xi_uc_S -1)
    Diff_c_K = 100*(xi_c_K / xi_uc_K -1)
    Diff_c_err = 0. # This is hard to calc. requires CxUC matrix; ignore for now. 

    # These first few ratios will
    # have UNCLIPPED FID on the denominator...:
    Diff_uc = np.zeros_like(xi_uc)
    Diff_S = np.zeros_like(xi_uc_S)
    Diff_err = np.zeros_like(xi_err)
    # (again, the sys, if dz, have different number of measurements per tomo bin,
    # so we have to have this special way of reading them in:)
    Diff_uc_Sys = []      
    Diff_uc_Sys_add = []  
    Diff_c_Sys_add = []   
    k = 0
    for i in range(len(noise)):
        for j in range(i, len(noise)):
            for c in range(len(cosmol)):
                Diff_uc[c,k,:] = 100*(xi_uc[c,k,:] / xi_uc[-1,k,:] -1)
            Diff_S[k,:] = 100*(xi_uc_S[k,:] / xi_uc[-1,k,:] -1)
            Diff_err[k,:] = 100*(xi_err[k,:] / xi_uc[-1,k,:])

            # dz X-bins have MORE bias measurements than auto-bins! SO we need to change array size:
            if Sys=="dz" and i!=j:
                files_array=Sys_files_all
            else:
                files_array=Sys_files

            Diff_uc_Sys_perz = np.zeros([len(files_array),nthbins]) 
            Diff_uc_Sys_add_perz = np.zeros([len(files_array),nthbins]) 
            Diff_c_Sys_add_perz = np.zeros([len(files_array),nthbins]) 

            for s in range(len(files_array)):
                Diff_uc_Sys_perz[s,:] = 100*(xi_uc_Sys2[k][s,:] / xi_uc[-1,k,:] -1)
                Diff_uc_Sys_add_perz[s,:] = 100*(xi_uc_Sys_add[k][s,:] / xi_uc[-1,k,:] -1)
                Diff_c_Sys_add_perz[s,:] = 100*(xi_c_Sys_add[k][s,:] / xi_c[-1,k,:] -1)
            Diff_uc_Sys.append( Diff_uc_Sys_perz )
            Diff_uc_Sys_add.append( Diff_uc_Sys_add_perz )
            Diff_c_Sys_add.append( Diff_c_Sys_add_perz )
            k += 1
    Diff_uc_K = xi_uc_K / xi_uc[-1]
    
    y_limits = [-129, 129]        
    # everything normalised to cosmoSLICS-fid:
    Ladder_Plot(avg_theta_S,
                Diff_uc,
                avg_theta_S,
                Diff_S,
                Diff_err,
                Diff_uc_Sys_add,
                #100*(Diff_uc_Sys - 1.),
                len(noise), r'$\theta \, [\rm{arcmin}]$',
                r'$\left( \xi_{%s} - \xi_{%s,{\rm fid}} \right) / \xi_{%s,{\rm fid}}$' %(pm,pm,pm) + ' [%]',
                [0.5, 300.], y_limits, None)


    y_limits_cuc = [[-999,999],
    [-99,99],
    [-99,99],
    [-99,99],
    [-99,99],
    [-99,99]]

    # clipped/unclipped:
    Ladder_Plot(avg_theta_S,
                Diff_c,
                avg_theta_S,
                Diff_c_S,
                Diff_c_err,
                Diff_c_Sys_add,
                len(noise), r'$\theta \, [\rm{arcmin}]$',
                r'$100 \times (\xi_{%s}^{\rm c} / \xi_{%s}^{\rm uc} -1)$' %(pm,pm),
                [0.5, 300.], y_limits_cuc, None)

else:
    # the sys doesn't transform easily, as an array, so:
    input_uc_sys = []
    input_c_sys = []
    for k in range(nzbins):
        input_uc_sys.append( avg_theta_S * 1e4 * xi_uc_Sys_add[k] )
        input_c_sys.append( avg_theta_S * 1e4 * xi_c_Sys_add[k] )

    # unclipped
    Ladder_Plot(avg_theta_S,
                avg_theta_S * 1e4 * xi_uc,
                avg_theta_S,
                avg_theta_S * 1e4 * xi_uc_S,
                avg_theta_S * 1e4 * xi_err,
                input_uc_sys,
                #avg_theta_S * 1e4 * xi_uc_Sys2,
                len(noise), r'$\theta \, [\rm{arcmin}]$',
                r'$\theta \times \xi_{%s}^{\rm uc} \, [10^{-4} \rm{arcmin}]$' %pm,
                [0.5, 300.], y_limits_uc, None)

    # clipped
    Ladder_Plot(avg_theta_S,
                avg_theta_S * 1e4 * xi_c,
                avg_theta_S,
                avg_theta_S * 1e4 * xi_c_S,
                avg_theta_S * 1e4 * xi_c_err,
                input_uc_sys,
                len(noise), r'$\theta \, [\rm{arcmin}]$',
                r'$\theta \times \xi_{%s}^{\rm c} \, [10^{-4} \rm{arcmin}]$' %pm,
                [0.5, 300.], y_limits_c, None)

# plot the raw predictions
#Ladder_Plot(avg_theta_CS,
#            xi_uc,
#            avg_theta_S,
#            xi_uc_S,
#            xi_err,
#            avg_theta_CS_fid_JHD,
#            xi_uc_fid_JHD,
#            len(noise), r'$\theta \, [\rm{arcmin}]$',
#            r'$\xi_{%s}$' %pm,
#            [0.5, 300.], ylog_limits, None)


# plot the ratio of the measurements to the fiducial cosmology
