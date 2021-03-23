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

- **processFilters.py**:

- **processSEDs.py**:

- **simulateWithSEDs.py**:



- **delightLearn.py**:	

							
- **templateFitting.py**:

- **delightApply.py**:							

## Other Delight modules

- **calibrateTemplateMixturePriors.py**: This module must be used outside RAI



## Python modules for interface purposes

-

- **makeConfigParam.py**:

- **convertDESCcat.py**:	

- **getDelightRedshiftEstimation.py**:

 







## Utilities

- **utils.py : **


