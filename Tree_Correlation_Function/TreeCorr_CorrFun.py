# 21/03/2017, B. M. Giblin, PhD student Edinburgh
# Create a config file for TreeCorr based on input parameters (KiDS or Mocks)
# And compute x+/- using TreeCorr

import treecorr
import time
import sys
import os.path
from os import getcwd
import numpy as np
from subprocess import call

pipeline_DIR='/home/bengib/Clipping_SimsLvW/'
data_DIR='/data/bengib/Clipping_SimsLvW/'
classdir = pipeline_DIR + "/ShowSumClass"

sys.path.insert(0, classdir) # add directory in which classes & functions 
							 # are defined to the python path

from ClassWarfare import Filter_Input


variable = Filter_Input(sys.argv)
variable.Filter()
RUN = sys.argv[1]

if RUN == 'Sims_Run':
	name, gpam, DIRname, SS, sigma, SN, mask, z, PS, sqdeg, zlo, zhi, ThBins, OATH, los, los_end = variable.Unpack_Sims()
	combined_name = '%s.%sGpAM.LOS%s'%(name,gpam,los)
	cn = los #config file designation
	if sqdeg == 60 or sqdeg == 100:
		metric='Euclidean'
		flip_g2=True
	elif sqdeg == 36:
		metric='Euclidean'
		flip_g2=False
	else:
		metric='Arc'	
		flip_g2=False # not certain if this is correct for Sims

else:
	DIRname, Blind, SS, sigma, zlo, zhi, ThBins, OATH, Field = variable.Unpack_KiDS()
	combined_name = '%s.Blind%s'%(Field,Blind)
	cn = Field # config file designation
	metric='Arc'
	flip_g2=True


if os.path.isdir('%s/Tree_Correlation_Function/%s/ThBins%s'%(data_DIR, DIRname,ThBins)) is False:
	call(['mkdir','-p', '%s/Tree_Correlation_Function/%s/ThBins%s'%(data_DIR, DIRname,ThBins)])



def Assemble_TreeCorr_ConfigFile(input_file, output_file, metric, flip_g2, bin_slop, min_sep, max_sep, ThBins, cn):
	f = open('%s/Tree_Correlation_Function/config_files/config_treecorr%s.yaml'%(pipeline_DIR,cn), 'w')
	f.write("#Configuration file for TreeCorr test" + "\n" + "#Input Parameters" + "\n" +
			"file_name: '%s'\n" %input_file +
			"ra_col: 1\n" + 
			"dec_col: 2\n" +
			"ra_units: 'degrees'\n" +
			"dec_units: 'degrees'\n" +
			"g1_col: 3\n" +
			"g2_col: 4\n" +
			"w_col: 5\n" +
			"metric: '%s'\n" %metric +
			"flip_g2: %s\n" %flip_g2 +
			"bin_slop: %s\n" %bin_slop +
			"#Output Parameters" + "\n" +
			"min_sep: %s\n" %min_sep +
			"max_sep: %s\n" %max_sep +
			"nbins: %s\n" %ThBins +
			"sep_units: 'arcmin'\n" +
			"gg_file_name: '%s'\n" %output_file) 
	f.close()
	return





min_sep=0.5
max_sep=300.
# Set bin_slop such that bin_size (in log space) * bin_slop <= 0.1 
# This is the minimum accuracy advised by Jarvis. Alternatively, set bin_slop = 0 for brute force
# OR bin_slop =< 0.1 which is indistinguishable from brute force in most cases (so says Jarvis in email to India).
bin_slop = 0.1/(np.log(max_sep/min_sep)/float(ThBins))


unclip_clip = [ 'ORIG', 'SS%s.rCLIP_%ssigma'%(SS,sigma) ]

for ucc in unclip_clip:

	# TreeCorr needs e2 flipped from what was correct for athena.
	# Made it so g2 is flipped in config_file
	input_file = '%s/Correlation_Function/%s/%s.%s.ThetaX_ThetaY_e1_e2_w.Std.asc' %(data_DIR, DIRname, combined_name, ucc)

	output_file = '%s/Tree_Correlation_Function/%s/ThBins%s/%s.%s.CorrFun.asc' %(data_DIR, DIRname, ThBins, combined_name, ucc)

	Assemble_TreeCorr_ConfigFile(input_file, output_file, metric, flip_g2, bin_slop, min_sep, max_sep, ThBins, cn)

	config_file='%s/Tree_Correlation_Function/config_files/config_treecorr%s.yaml' %(pipeline_DIR,cn)
	config = treecorr.read_config(config_file)

	t1 = time.time()
	treecorr.corr2(config)
	t2=time.time()

	print("TreeCorr time for %s is %.1f s" %(output_file, (t2-t1)))



