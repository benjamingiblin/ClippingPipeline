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



## The most important codes at a glance

** (0) Launch.sh **

The slurm script which launches the clipping pipeline on cuillin. This asks for the clipping pipeline job to be put in the queue to run on one of the cuillin workers. I have set the memory request of the launch script low enough so that the clipping pipeline will get to the front of the queue quickly, but not so low that the job fails. 

** (1) Master_CorrFun_ByParts.sh **

The master bash script that runs the whole pipeline. It takes an input parameter file which specifies the details of the calculation (more on this below). The master script first of all executes scripts which find the relevant input shear catalogues specified by the parameter file, sets up temporary directories on the cuillin workers to save outputs to, and copies some files over to the cuillin worker. It then starts running the pipeline.

** (2) SkelePipeline_PartI.sh and SkelePipeline_PartII.sh **
The two bash scripts, executed by Master_CorrFun_ByParts.sh, which run the first part and second parts of the pipeline respectively.

SkelePipeline_PartI.sh runs the following:
 * KiDSMocks_DataGrab.sh/KiDS450_DataGrab.sh - accesses the input simulated/data shear catalogues.
 * Sims_DataGrab.py - performs redshift cuts, injects random galaxy shape noise to the simulated catalogues.
 * Mass_Recon/MassMapLvW.sh - performs the shear-to-kappa mass reconstruction.
 * Clipping_K/get_clipped_shear.py - clips the kappa maps, converts kappa back to shear, and saves a clipped shear catalogue.

SkelePipeline_PartII.sh runs the following:
 * Correlation_Function/Convert_2_AthenaStandard.sh - converts the clipped shear catalogue, and the original unclipped shear catalogues into a format which is understood by the code to calculate the shear correlation function (TreeCorr).
 * Tree_Correlation_Function/TreeCorr_CorrFun.py - calculates the clipped and unclipped shear correlation functions and saves the output.


## More detail on the codes

** (4) 