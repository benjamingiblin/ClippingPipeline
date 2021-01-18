# 18/01/2020, B. M. Giblin, Postdoc, Uni. Edinburgh

# The Clipping Pipeline

The following is a README document for the clipping pipeline employed in Giblin et al. (2020) - https://arxiv.org/abs/1805.12084

This series of codes start from shear catalogues perform the following, either with simulations or real data:

 * Beginning with shear catalogues, redshift cuts are made, and intrinsic galaxy shape noise is optionally added if working with simulations.
 * The shear catalogues are converted into 2D convergence (projected surface density) maps using the flat-sky Kaiser-Squires '93 calculation. The maps are smoothed with a Gaussian filter with a designated standard deviation (smoothing scale).
 * The convergence maps are "clipped" at a designated threshold: any pixels where the convergence exceeds the threshold, the values get set equal to the threshold.
 * The clipped convergence map is subtracted from the original "unclipped" map. This creates the "residual" convergence map, containing only the peaks above the threshold, and zeroes everywhere else.
 * We convert the residual convergence map back to "residual" shear by performing the inverse of the Kaiser-Squires calculation. This generates a map of the residual shear signal corresponding to the peaks. We convert this residual shear map into a catalogue.
 * The residual shear is subtracted from the original ("unclipped") shear, creating the "clipped" shear catalogue.
 * From this catalogue, we calculate the clipped shear correlation function using TreeCorr (credit to Mike Jarvis).
 * The unclipped shear catalogue is also fed to TreeCorr to calculate the conventional, "unclipped" shear correlation functions.


## General line used to execute the pipeline:

For simulations
./Master_CorrFun_ByParts.sh Sims_Run param_files/<input parameter file> <los_start> <los_end>
where one simply needs to designate the input parameter file, and the line of sight numbers (IDs of the simulation catalogues) to start and end the clipping on.

Note that on the supercomputer cuillin,	you will be having cuillin run the pipeline on its various workhorse computers (hereafter, referred to simply as "workers"). This means you will be changing that line above in the Launch.sh script (explained in the next section), and running it repeatedly.


## Installation

For this pipeline to work, you need to have already installed an anaconda python distribution from, e.g., https://docs.anaconda.com/anaconda/install/

Then, follow these steps to install:





## The most important codes at a glance

*(0) Launch.sh*

 * Executed as: sbatch Launch.sh 

This script "launches" the clipping pipeline on cuillin. This means that the script asks for the clipping pipeline job to be put in the queue to run on one of the cuillin workers. I have set the memory request of the launch script low enough so that the clipping pipeline will get to the front of the queue quickly, but not so low that the job fails. 

*(1) Master_CorrFun_ByParts.sh*

 * Executed as ./Master_CorrFun_ByParts.sh Sims_Run param_files/<input parameter file> <los_start> <los_end>
 * (note this will only work on a cuillin worker, not on the head node itself).

The master bash script that runs the whole pipeline. It takes an input parameter file which specifies the details of the calculation (more on this below). The master script first of all executes scripts which find the relevant input shear catalogues specified by the parameter file, sets up temporary directories on the cuillin workers to save outputs to, and copies some files over to the cuillin worker. It then starts running the pipeline.

*(2) SkelePipeline_PartI.sh and SkelePipeline_PartII.sh*
The two bash scripts, executed by Master_CorrFun_ByParts.sh, which run the first part and second parts of the pipeline respectively.

 * Executed as ./SkelePipeline_PartI.sh Sims_Run param_files/<input parameter file> <los_start> <los_end> 
    
SkelePipeline_PartI.sh runs the following:
 * KiDSMocks_DataGrab.sh/KiDS450_DataGrab.sh - accesses the input simulated/data shear catalogues and applies redshift cuts.
 * Sims_DataGrab.py - injects random galaxy shape noise to the simulated catalogues.
 * Mass_Recon/MassMapLvW.sh - performs the shear-to-kappa mass reconstruction.
 * Clipping_K/get_clipped_shear.py - clips the kappa maps, converts kappa back to shear, and saves a clipped shear catalogue.

SkelePipeline_PartII.sh runs the following:
 * Correlation_Function/Convert_2_AthenaStandard.sh - converts the clipped shear catalogue, and the original unclipped shear catalogues into a format which is understood by the code to calculate the shear correlation function (TreeCorr).
 * Tree_Correlation_Function/TreeCorr_CorrFun.py - calculates the clipped and unclipped shear correlation functions and saves the output.


## More detail on the codes

*(4) KiDSMocks_DataGrab.sh*

This code pulls the (RA, Dec, shear1, shear2, weight, z_spec, z_phot) information out of the whatever shear catalogue is specified by the input parameter file. This code applies the specified redshift cuts and saves an output shear catalogue. 



