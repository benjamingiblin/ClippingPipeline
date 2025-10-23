# 21/08/2025
# make nice PDF plot for ERC app

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

# just use highest signal bin
SN = '0.254'
ZB = '0.7-0.9'
DIR = 'MRres140.64arcs_100Sqdeg_SN%s_Mosaic_KiDS1000GpAM_zKiDS1000_ZBcut%s_Cosmol' %(SN,ZB)

nkbins = 20
SS = 2.816
cosmol = []
for i in range(25):
    cosmol.append('%s' %i)
cosmol.append('fid')

PDFs = np.zeros([ len(cosmol), nkbins ])

for i in range(len(cosmol)):
    filename = '%s%s/SN%s_Mosaic.KiDS1000GpAM.LOSAll.SS%s.SNRPDF_%sbins.dat' %(DIR,cosmol[i],SN,SS,nkbins)
    snr, PDFs[i,:] = np.loadtxt(filename, usecols=(0,1),unpack=True)

DIR_slics = DIR.replace('_Cosmol','').replace('SN%s'%SN,'SNCycle')
filename_slics = '%s/SN%s_Mosaic.KiDS1000GpAM.NLOS2170-Shuffle0.SS%s.SNRPDF_%sbins.CovMat.npy' %(DIR_slics,SN,SS,nkbins)
err = np.sqrt( np.diag( np.load(filename_slics) ))
    
filename_kids = filename.replace('Cosmolfid','CosmolKiDS1000').replace('SN%s'%SN,'SNKiDS')
PDF_kids = np.loadtxt(filename_kids, usecols=(1,), unpack=True)


S8 = np.loadtxt('/home/bengib/cosmoSLICS/Cosmologies/SLICS_Cosmol_Table_ExtraCol_25plusFid.txt',
                usecols=(1,), unpack=True)
norm = matplotlib.colors.Normalize(vmin=S8.min(), vmax=S8.max())
c_m = matplotlib.cm.jet #cool
s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
s_m.set_array([])

idx = np.where(abs(snr)<2.3)[0] # just plot central SNR's
fig = plt.figure(figsize = (10,8.5))
for i in range(len(cosmol)):
    plt.plot(snr[idx], (PDFs[i]/PDFs[-1])[idx], color=s_m.to_rgba(S8[i]), linewidth=2 )
plt.errorbar(snr[idx], (PDF_kids/PDFs[-1])[idx], (err/PDFs[-1])[idx],
             color='magenta', linewidth=3)
plt.xlabel('Density/Noise')
plt.title('PDF / PDF$(S_8=%.2f)$'%S8[-1])
plt.ylim([0.88, 1.07])
plt.xlim([-2.4,2.4])
fig.subplots_adjust(hspace=0, wspace=0)
fig.subplots_adjust(right=0.7)
cbar_ax = fig.add_axes([0.72, 0.15, 0.05, 0.7])
fig.colorbar(s_m, cax=cbar_ax, label=r'$S_8=\sigma_8 \sqrt{\Omega_{\rm m} / 0.3}$' )
plt.savefig('Figures_4_Paper/ERC_PDFs.png')
plt.show()

    
