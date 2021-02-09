# 18/01/2020, B. M. Giblin, Postdoc, Uni. Edinburgh

# The Clipping Pipeline

The following is a README document for the clipping pipeline employed in Giblin et al. (2018) - https://arxiv.org/abs/1805.12084

This series of codes start from shear catalogues perform the following, either with simulations or real data:

 * Beginning with shear catalogues, redshift cuts are made, and intrinsic galaxy shape noise is optionally added if working with simulations.
 * The shear catalogues are converted into 2D convergence (projected surface density) maps using the flat-sky Kaiser-Squires '93 calculation. The maps are smoothed with a Gaussian filter with a designated standard deviation (smoothing scale).
 * The convergence maps are "clipped" at a designated threshold: any pixels where the convergence exceeds the threshold, the values get set equal to the threshold.
 * The clipped convergence map is subtracted from the original "unclipped" map. This creates the "residual" convergence map, containing only the peaks above the threshold, and zeroes everywhere else.
 * We convert the residual convergence map back to "residual" shear by performing the inverse of the Kaiser-Squires calculation. This generates a map of the residual shear signal corresponding to the peaks. We convert this residual shear map into a catalogue.
 * The residual shear is subtracted from the original ("unclipped") shear, creating the "clipped" shear catalogue.
 * From this catalogue, we calculate the clipped shear correlation function using TreeCorr (credit to Mike Jarvis).
 * The unclipped shear catalogue is also fed to TreeCorr to calculate the conventional, "unclipped" shear correlation functions.



## Installation

For this pipeline to work, you need to have already installed an anaconda python distribution from, e.g., https://docs.anaconda.com/anaconda/install/

You also need to have the treecorr package (https://github.com/rmjarvis/TreeCorr). Once anaconda is installed, TreeCorr should be easy to install too, with: pip install treecorr

Then, follow these steps to install the clipping pipeline. For the purposes of illustration, these instructions describe how to install in your home directory on cuillin, although the pipeline could happily be installed elsewhere on this machine.

1. Navigate to your home directory ( /home/<user_name> ), then git clone the pipeline:

   git clone https://github.com/benjamingiblin/ClippingPipeline.git Clipping_Pipeline

2. Move the Install_Pipeline subdirectory upwards into the home directory:

   mv Clipping_Pipeline/Install_Pipeline .

3. Open Install_Pipeline/Install_Pipeline.sh. Comment-out the fail-safe line "exit 0" which prevents the installation script being ran accidentally. Ensure that the Pipeline_DIR variable correctly points to the subdirectory into which you cloned the pipeline. Assuming you named the subdirectory "Clipping_Pipeline" as suggested by step 1., then this line shouldn't need any modification.

4. Run the installation script:

   Install_Pipeline/Install_Pipeline.sh

5. Open your .bashrc script (on cuillin, this is stored in /home/<user_name>/). Copy and paste the following line at the bottom: ulimit -s unlimited


Note that your .bashrc should be automatically sourced every time you log on to cuillin. If it isn't, then anaconda will not be set to your default python package and your ulimit will not be set to unlimited every time you log on. If you're .bashrc is not being automatically sourced (check your python path to confirm if it is), then you can correct this by making a .bash_profile script in your home diretcory containing this line only:

source ~/.bashrc

Also note it is recommended you remove the commenting on the "exit 0" line in the installation script after this is done, to prevent it being ran accidentally.

The pipeline should now be correctly configured for the user. 




## Executing the pipeline

The command to run the pipeline is:

./Master_CorrFun_ByParts.sh Sims_Run param_files/<input_parameter_file> <los_start> <los_end>

where one simply needs to designate the input parameter file, and the line of sight numbers (IDs of the simulation catalogues) to start and end the clipping on.


Note that the pipeline is designed to run on the workhorse processors (hereafter referred to as "workers") of the cuillin supercomputer at the ROE. This means that the pipeline will not run successfully if executed on the cuillin head node. Running the pipeline needs to be designated to the workers using the launch script via:

sbatch Launch_Pipeline.sh

In order to run different pipeline settings, one simply needs to change the <input_parameter_file>, <los_start> and <los_end> variables in the Launch_Pipeline.sh script.

If you want to quickly check the pipeline runs successfully without launching a job, waiting for it to be allocated to a worker and finishing running, you can execute directly on the cuillin worker manually. To do this, ssh into a worker (with, e.g., ssh worker019), navigate into the Clipping_Pipeline directory and run the line to execute the pipeline there. NOTE: you should only do this for 1-3 lines of sight in total, which will only take a few minutes to run. Taking up processing power and memory  on a worker for long durations without using the launch script is in general bad practice.



## The input parameter file

These are stored in the param_files/ subdirectory. Here is the contents of a parameter file for clipping the cosmoSLICS simulations which have been tailored to match the KiDS1000 data set. The most important parameters which the user may want to change are commented on below.

     KiDS1000			# number of gals/arcmin^2 (KiDS1000 means set equal to the KiDS1000 data)
     9.33 			# Smoothing scale [pxls] x sqrt(2) 
     X3 			# The clip threshold (see Clipping_K/Clip_Thresholds)
     0.27 			# The intrinsic ellipticity of galaxies, to be added if running on simulations.
     nomask			# Whether masking is to be included or omitted.
     KiDS1000 			# The n(z) of the galaxies in the data (KiDS1000 means set to the SOM-GOLD KiS1000 n(z) )
     0.00129115558		# Pixel scale of the mocks (deg/pxl0 - REDUNDANT.
     100			# Angular size of the mocks [deg^2]
     None			# Cosmology ID. For cosmoSLICS [fid,0-24]. For SLICS, put None.
     0.1			# zlow: lower limit on zB cut. For no zB cut, put anything that isn't a number
     0.3			# zhigh: upper limit on zB cut. For no zB cut, put anything that isn't a number	
     9				# The number of theta bins to calculate the xi_+/- in [uses log-spaced bins 0.5-300arcmin]
     60arcs		        # Resolution of the kappa maps: 60arcs = 60 arcsecons per pxl. 
     0.02			# OATH value previously used in Athena [Redundant]. 


**The most important parameters**

	- ROW 2: 9.33 is the smoothing scale (SS) x sqrt(2) and given in units of pixels. The SS controls how much the maps are smoothed in mass reconstruction (specifically it is the standard deviation of the Gaussian smoothing filter). To convert the value in the param-file to an angular smoothing scale, use the resolution of the map set in the penultimate row. Here the resolution is 60 arcseconds/pxl, so SS*sqrt(2)=9.33 corresponds to 396 arcsec = 6.6 arcmin.
	- ROW 3: X3 designates the clipping threshold. This tells the pipeline to read the convergence value saved in the file: Clipping_K/Clip_Thresholds/X3_Threshold (which is 0.01). Other thresholds (X0-X4) can be found in this subdirectory, and the user is free to create new threshold files and use those.
	- ROW 4: 0.27 is the standard deviation of the Gaussian shape noise (SN) added to the shear values in simulations. For no shape noise, this should be set to 0.
	- ROWS 10 AND 11: 0.1 and 0.3: These are the lower and upper redshift cuts to make on the data. Working with KiDS1000-like simulations, these numbers should always be (0.1-0.3, 0.3-0.5, 0.5-0.7, 0.7-0.9, 0.9-1.2). Note that the SN value is different for these 5 redshift bins: (0.27, 0.258, 0.273, 0.254, 0.27)





## The most important codes at a glance

*(0) Launch_Pipeline.sh*

 * Executed as: sbatch Launch_Pipeline.sh 

This script "launches" the clipping pipeline on cuillin. This means that the script asks for the clipping pipeline job to be put in the queue to run on one of the cuillin workers. I have set the memory request of the launch script low enough so that the clipping pipeline will get to the front of the queue quickly, but not so low that the job fails. 

*(1) Master_CorrFun_ByParts.sh*

 * Executed as ./Master_CorrFun_ByParts.sh Sims_Run param_files/<input_parameter_file> <los_start> <los_end>
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



