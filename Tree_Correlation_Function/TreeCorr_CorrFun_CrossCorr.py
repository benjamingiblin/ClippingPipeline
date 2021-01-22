# 21/03/2017, B. M. Giblin, PhD student Edinburgh
# Compute the cross-correlation between different redshifts
# INPUT IS: Sims/KiDS_RUN paramfile1 los los DUMMY Sims/KiDS_RUN paramfile2 los los
# DUMMY is important keyword that muct be included....^


import treecorr
import time
import sys
import os.path
from os import getcwd
import numpy as np
from subprocess import call

overall_DIR = getcwd()
classdir = getcwd() + "/ShowSumClass"
sys.path.insert(0, classdir) # add directory in which classes & functions 
							 # are defined to the python path

from ClassWarfare import Filter_Input

print( "Reading and filtering inputs for first file" )
#print sys.argv[0:5]
variable = Filter_Input(sys.argv[0:5])
variable.Filter()
RUN = sys.argv[1]

print( "Reading and filtering inputs for second file" )
#print sys.argv[5:]
variable2 = Filter_Input(sys.argv[5:])
variable2.Filter()




if RUN == 'Sims_Run':
        name1, gpam, DIRname1, SS, sigma, SN, mask, z, PS, sqdeg, zlo1, zhi1, ThBins, OATH, los, los_end = variable.Unpack_Sims()
        name2, gpam, DIRname2, SS, sigma, SN, mask, z, PS, sqdeg, zlo2, zhi2, ThBins, OATH, los, los_end = variable2.Unpack_Sims()

        combined_name1 = '%s.%sGpAM.LOS%s'%(name1,gpam,los)
        combined_name2 = '%s.%sGpAM.LOS%s'%(name2,gpam,los)
        cn = los #config file designation
        if sqdeg == 60 or sqdeg == 100:
                metric='Euclidean'
        else:
                metric='Arc'	
        flip_g2=True # not certain if this is correct for Sims

else:
        DIRname1, Blind, SS, sigma, zlo1, zhi1, ThBins, OATH, Field = variable.Unpack_KiDS()
        DIRname2, Blind, SS, sigma, zlo2, zhi2, ThBins, OATH, Field = variable2.Unpack_KiDS()
        combined_name = '%s.Blind%s'%(Field,Blind)
        cn = Field # config file designation
        metric='Arc'
        flip_g2=True



# Make the X-correlation output directory
ZBcut2 = DIRname2.split('ZBcut')[-1]
DIRname = '%s_X_ZBcut%s' %(DIRname1, ZBcut2)
if os.path.isdir('%s/Tree_Correlation_Function/%s/ThBins%s'%(overall_DIR, DIRname,ThBins)) is False:
	call(['mkdir','-p', '%s/Tree_Correlation_Function/%s/ThBins%s'%(overall_DIR, DIRname,ThBins)])



def Assemble_TreeCorr_ConfigFile(input_file1, input_file2, output_file, metric, flip_g2, bin_slop, min_sep, max_sep, ThBins, cn):
	f = open('%s/Tree_Correlation_Function/config_files/config_treecorr%s.yaml'%(overall_DIR,cn), 'w')
	f.write("#Configuration file for TreeCorr test" + "\n" + "#Input Parameters" + "\n" +
			"file_name: [ %s, %s ]\n" %(input_file1, input_file2) +
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


#for ucc in unclip_clip:
for cycle in range(3):

        if cycle == 0:
                # Do unclipped X unclipped
                keyword1 = 'ORIG'
                keyword2 = keyword1
                out_keyword  = keyword1
                
        elif cycle == 1:
                # Do clipped X clipped
                keyword1 = 'SS%s.rCLIP_%ssigma'%(SS,sigma)
                keyword2 = keyword1
                out_keyword  = keyword1

        elif cycle == 2:
                # Do unclipped X clipped
                keyword1 = 'ORIG'
                keyword2 = 'SS%s.rCLIP_%ssigma'%(SS,sigma)
                out_keyword  = keyword1 + '_X_' + keyword2

        input_file1 = '%s/Correlation_Function/%s/%s.%s.ThetaX_ThetaY_e1_e2_w.Std.asc' %(overall_DIR, DIRname1, combined_name1, keyword1)
        input_file2 = '%s/Correlation_Function/%s/%s.%s.ThetaX_ThetaY_e1_e2_w.Std.asc' %(overall_DIR, DIRname2, combined_name2, keyword2)	
        output_file = '%s/Tree_Correlation_Function/%s/ThBins%s/%s.%s.CorrFun.asc' %(overall_DIR, DIRname, ThBins, combined_name1, out_keyword)

        Assemble_TreeCorr_ConfigFile(input_file1, input_file2, output_file, metric, flip_g2, bin_slop, min_sep, max_sep, ThBins, cn)

        config_file='%s/Tree_Correlation_Function/config_files/config_treecorr%s.yaml' %(overall_DIR,cn)
        config = treecorr.read_config(config_file)

        t1 = time.time()
        treecorr.corr2(config)
        t2=time.time()

        print( "TreeCorr time for %s is %.1f s" %(output_file, (t2-t1)) )

               
