# README.md : Interface of DESC-LSST RAIL to Delight package 

---------------------------------------------------------------------------------------

## RAIL interface class to Delight

Steering Delight from RAIL


Used on Vera C. Runbin LSST only estimation

- Author        : Sylvie Dagoret-Campagne 
- Affiliation   : IJCLab/IN2P3/CNRS/France
- Creation date : March 23th 2021
- Last update   : March 24th 2021

----------------------------------------------------------------------------------------



These set of python modules may be called from RAIL **DelightPZ** algorithm for PhotoZ estimation with Delight.
When Delight is used in the context of DESC-LSST, a training dataset with a validation dataset of photometrics redshift data are provided in RAIL.
When running Delight in RAIL-DESC context, we are refering to the running **STANDARD MODE**.

Optionaly it is possible to run Delight in **TUTORIAL MODE**. In the later case Delight generate internally training and test datasets. However, there is no connexion data between RAIL and Delight.
One has to consider the **TUTORIAL MODE** should be used by expert only for internal control and debugging purposes.

Delight is controled by its configuration parameter file. This means that this configuration parameter file is generated from RAIL according configuration parameters set in delightPZ.yaml file.

Then the RAIL traning and validation dataset in DC2 (hpf5) are converted into an ascii format file readable by Delight. Conversely results of PhotoZ estimation are written in an ascii files, then must be given to RAIL. 

## Main Delight scripts transposed to be called by rail

### Initialisation

The two following modules are necessary in **STANDARD MODE**. They requires the Filter dataset and SED datasets provided by Delight.
By now, the Delight setpy.py has been modified to ship those data in the installation directory of Delight.

- **processFilters.py**: Decompose the Filter transmission into a gaussian mixture sum. The filters to be used are selected in the delight configuration file. In the context of DESC, only LSST filters are selected (u,g,r,i,z,y). The parameters of the fit are written in the file '_gaussian_coefficients.txt' in the Filter direcptry.

- **processSEDs.py**: Build a model of fluxes-redshifts from a set of SED templates in the chosen filters. It is possible to use either the CWW model or the Brown model. By now, CWW SED are chosen in the delight configuration file. The fluxes-redshifts models, one per SED template are written in the file '_fluxredshiftmod.txt'.

### creation

- **simulateWithSEDs.py**: In TUTORIAL MODE only, it generates mock dataset inside Delight. These data are the traning and validation samples of fluxes-fluxeserrors-redshift drawing ramdomly the type of the underlying SED template Fluxes-Redshifts model build in **processSEDs.py**. These data are written in the ascii files 'trainingFile' and 'targetFile' which true name is given in the configuration parameter file.

### Learn

- **delightLearn.py**:	Build the Gaussian Process from the training data. The Gaussian process is written in the file 'training_paramFile'.

### Estimation	

The next modules generate the pdf of the redshift estimate for each sample of the validation dataset (not filtered out ones).
							
- **templateFitting.py**: Does the estimation with template fitting. It writes the file with the z posterior in 'redshiftpdfFileTemp'.

- **delightApply.py**: Does the estimation with template fitting. It writes the file with the z posterior in write the 'redshiftpdfFile'.			

## Other Delight modules

- **calibrateTemplateMixturePriors.py**: This module must be used outside RAIL.


## Python modules for interface purposes


- **makeConfigParam.py**: Generate Delight configuration file from Delight configuration file. In the case of RAIL estimation, the validation data is split into different chunk. Then there must be is one Delight configuration file per chunk. A Delight configuration file must contain all data file paths used by Delight.

- **convertDESCcat.py**: Convert DESC PhotoZ dataset into ascii file readable by Delight. For now, hdf5 file format are used. It is important to notice that for DESC, some input data filering has been applied. For example some of the samples having bad fluxes, or missing fluxes bands are removed from training and validation dataset. Moreover these data can be filtered with bad fluxes SNR.
This selection is controled from the RAIL delightPZ.yaml file.
This module handle the validation dataset split into different chunks.

- **getDelightRedshiftEstimation.py**: Retrieve results on PhotoZ estimation of Delight for RAIL. For each chunk of the validation dataset, the mode of redshift distribution and its error are returned.
For those entries  that have been filtered by **convertDESCcat.py**, a -1 value is returned(meaning bad sample).





## Utilities

- **utils.py** : Library of IO inside RAIL used in **convertDESCcat.py** for reading DESC dataset.


## Status of the code

The Delight / RAIL interfacing code is quite in progress. The biggest work remaining to be done is deep Delight optimization with DESC-DC2.



## Example of Delight parameter file

Below are the configuration file parameters required by Delight. This parameter file is generated by **makeConfigParam.py**.
Some parameters are dedined in delightPZ.yaml, some others are hardcoded.
Later, some of the hardcoded parameters will be moved to be tuned in delightPZ.yaml.


	[Bands]
	names: lsst_u lsst_g lsst_r lsst_i lsst_z lsst_y
	directory: /home/ubuntu/.local/lib/python3.8/site-packages/delight/../data/FILTERS
	numCoefs: 15
	bands_verbose: True
	bands_debug: True
	bands_makeplots: False

	[Templates]
	directory: /home/ubuntu/.local/lib/python3.8/site-packages/delight/../data/CWW_SEDs
	names: El_B2004a Sbc_B2004a Scd_B2004a SB3_B2004a SB2_B2004a Im_B2004a ssp_25Myr_z008 ssp_5Myr_z008
	p_t: 0.27 0.26 0.25 0.069 0.021 0.11 0.0061 0.0079
	p_z_t:0.23 0.39 0.33 0.31 1.1 0.34 1.2 0.14
	lambdaRef: 4.5e3

	[Simulation]
	numObjects: 1000
	noiseLevel: 0.03
	trainingFile: ./tmp/delight_data/galaxies-fluxredshifts.txt
	targetFile: ./tmp/delight_data/galaxies-fluxredshifts2.txt
	
	[Training]
	catFile: ./tmp/delight_data/galaxies-fluxredshifts.txt
	bandOrder: lsst_u lsst_u_var lsst_g lsst_g_var lsst_r lsst_r_var lsst_i lsst_i_var lsst_z lsst_z_var lsst_y lsst_y_var redshift
	referenceBand: lsst_i
	extraFracFluxError: 1e-4
	paramFile: ./tmp/delight_data/galaxies-gpparams.txt
	crossValidate: False
	CVfile: ./tmp/delight_data/galaxies-gpCV.txt
	crossValidationBandOrder: _ _ _ _ lsst_r lsst_r_var _ _ _ _ _ _
	numChunks: 1

	[Target]
	catFile: ./tmp/delight_data/galaxies-fluxredshifts2.txt
	bandOrder: lsst_u lsst_u_var lsst_g lsst_g_var lsst_r lsst_r_var lsst_i lsst_i_var lsst_z lsst_z_var lsst_y lsst_y_var redshift
	referenceBand: lsst_r
	extraFracFluxError: 1e-4
	redshiftpdfFile: ./tmp/delight_data/galaxies-redshiftpdfs.txt
	redshiftpdfFileTemp: ./tmp/delight_data/galaxies-redshiftpdfs-cww.txt
	metricsFile: ./tmp/delight_data/galaxies-redshiftmetrics.txt
	metricsFileTemp: ./tmp/delight_data/galaxies-redshiftmetrics-cww.txt
	useCompression: False
	Ncompress: 10
	compressIndicesFile: ./tmp/delight_data/galaxies-compressionIndices.txt
	compressMargLikFile: ./tmp/delight_data/galaxies-compressionMargLikes.txt
	redshiftpdfFileComp: ./tmp/delight_data/galaxies-redshiftpdfs-comp.txt
	
	[Other]
	rootDir: ./
	zPriorSigma: 0.2
	ellPriorSigma: 0.5
	fluxLuminosityNorm: 1.0
	alpha_C: 1.0e3
	V_C: 0.1
	alpha_L: 1.0e2
	V_L: 0.1
	lines_pos: 6500 5002.26 3732.22
	lines_width: 20.0 20.0 20.0
	redshiftMin: 0.01
	redshiftMax: 3.01
	redshiftNumBinsGPpred: 300
	redshiftBinSize: 0.001
	redshiftDisBinSize: 0.2
	confidenceLevels: 0.1 0.50 0.68 0.95





