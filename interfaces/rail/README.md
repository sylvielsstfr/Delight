# README.md : Interface of DESC-LSST RAIL to Delight package 

---------------------------------------------------------------------------------------

## RAIL interface class to Delight

Steering Delight from RAIL


Used on Vera C. Runbin LSST only estimation

- Author        : Sylvie Dagoret-Campagne

- Affiliation   : IJCLab/IN2P3/CNRS/France

- Creation date : March 2021

- Last update   : March 23th 2021

----------------------------------------------------------------------------------------



These set of python modules may be called from RAIL **DelightPZ** algorithm for PhotoZ estimation with Delight.
When Delight is used in the context of DESC-LSST, a training dataset with a validation dataset of photometrics redshift data are provided in RAIL.
This is called the **STANDARD MODE**.

Optionaly it is possible to run Delight in **TUTORIAL MODE**. In the later case Delight generate internally training and test datasets. However, there is no connexion data between RAIL and Delight.
One has to consider the **TUTORIAL MODE** should be used by expert only for debugging purposes.

Delight is controled by its configuration parameter file. This means that this configuration parameter file is generated from RAIL according configuration parameter.

Then the RAIL traning and validation dataset are converted into an ascii format file readable by Delight. Conversely result of PhotoZ estimation are written in an ascii file, then must be given to TAI. 

## Main Delight scripts transposed to be called by rail

### Initialisation

- **processFilters.py**: Decompose the Filter transmission into a gaussian mixture sum. The filters to be used are selected in the delight configuration file. In the context of DESC, only LSST filters are selected (u,g,r,i,z,y). The parameters of the fit are written in the file '_gaussian_coefficients.txt' in the Filter direcptry.

- **processSEDs.py**: Build a model of fluxes-redshifts from a set of SED templates in the chosen filters. It is possible to use either the CWW model or the Brown model. By now, CWW SED are chosen in the delight configuration file. The fluxes-redshifts models, one per SED template are written in the file '_fluxredshiftmod.txt'.

### creation

- **simulateWithSEDs.py**: In TUTORIAL MODE generate mock dataset inside Delight. These data are the traning and validation samples of fluxes-fluxeserrors-redshift drawing ramdomly the type of the underlying SED template Fluxes-Redshifts model build in **processSEDs.py**. These data are written in the ascii files 'trainingFile' and 'targetFile' which true name is given in the configuration parameter file.

### Learn

- **delightLearn.py**:	Build the Gaussian Process from the training data. The Gaussian process is written in the file 'training_paramFile'.

### Estimation	
							
- **templateFitting.py**:

- **delightApply.py**:							

## Other Delight modules

- **calibrateTemplateMixturePriors.py**: This module must be used outside RAIL.


## Python modules for interface purposes


- **makeConfigParam.py**: Generate Delight configuration file from Delight configuration file. In the case of RAIL estimation, the validation data is split into different chunk. Then there must be is one Delight configuration file per chunk. A Delight configuration file must contain all data file paths used by Delight.

- **convertDESCcat.py**: Convert DESC PhotoZ dataset into ascii file readable by Delight. For now, hdf5 file format are used.

- **getDelightRedshiftEstimation.py**: Retrieve results on PhotoZ estimation of Delight for RAIL. 

 



## Utilities

- **utils.py** : Library of IO inside RAIL used in **convertDESCcat.py** for reading DESC dataset


