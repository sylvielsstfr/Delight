#######################################################################################################
#
# script : convertDESCcat.py
#
# convert DESC catalog to be injected in Delight
# produce files `galaxies-redshiftpdfs.txt` and `galaxies-redshiftpdfs2.txt` for training and target
#
#########################################################################################################


import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from delight.io import *
from delight.utils import *
import h5py

logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG', logger=logger,fmt='%(asctime)s,%(msecs)03d %(programname)s, %(name)s[%(process)d] %(levelname)s %(message)s')


def convertDESCcat(configfilename,desctraincatalogfile,desctargetcatalogfile):
    """

    Convert files in ascii format to be used by Delight

    input args:
    - configfilename : Delight configuration file containg path for output files (flux variances and redshifts)
    - desctraincatalogfile : training file provided by RAIL (hdf5 format)
    - desctargetcatalogfile : target file provided by RAIL (hdf5 format)

    output :
    the Delight training and target file which path is in configuration file

    :param configfilename:
    :return:
    """


    logger.info("--- Convert DESC training and target catalogs ---")

    # 1) DESC catalog file
    msg="read DESC hdf5 file {} ".format(desctraincatalogfile)
    logger.debug(msg)

    with h5py.File(desctraincatalogfile, "r") as f:
        # List all groups
        msg="Keys: %s" % f.keys()
        logger.debug(msg)
        a_group_key = list(f.keys())[0]

        # Get the data
        descdata = list(f[a_group_key])

    logger.debug(descdata)

    f = h5py.File(desctraincatalogfile, 'r')

    galid = f['photometry/id'][()]
    redshifts = f['photometry/redshift'][()]

    mag_err_g_lsst = f['photometry/mag_err_g_lsst'][()]
    mag_err_i_lsst = f['photometry/mag_err_i_lsst'][()]
    mag_err_r_lsst = f['photometry/mag_err_r_lsst'][()]
    mag_err_u_lsst = f['photometry/mag_err_u_lsst'][()]
    mag_err_y_lsst = f['photometry/mag_err_y_lsst'][()]
    mag_err_z_lsst = f['photometry/mag_err_z_lsst'][()]
    mag_g_lsst = f['photometry/mag_g_lsst'][()]
    mag_i_lsst = f['photometry/mag_i_lsst'][()]
    mag_r_lsst = f['photometry/mag_r_lsst'][()]
    mag_u_lsst = f['photometry/mag_u_lsst'][()]
    mag_y_lsst = f['photometry/mag_y_lsst'][()]
    mag_z_lsst = f['photometry/mag_z_lsst'][()]

    # conversion of magnitude to flux
    flux_g_lsst = np.power(10, -0.4 * mag_g_lsst)
    flux_i_lsst = np.power(10, -0.4 * mag_i_lsst)
    flux_r_lsst = np.power(10, -0.4 * mag_r_lsst)
    flux_u_lsst = np.power(10, -0.4 * mag_u_lsst)
    flux_y_lsst = np.power(10, -0.4 * mag_y_lsst)
    flux_z_lsst = np.power(10, -0.4 * mag_z_lsst)

    # Flux error
    flux_err_u = mag_err_u_lsst * flux_u_lsst
    flux_err_g = mag_err_g_lsst * flux_g_lsst
    flux_err_r = mag_err_r_lsst * flux_r_lsst
    flux_err_i = mag_err_i_lsst * flux_i_lsst
    flux_err_z = mag_err_z_lsst * flux_z_lsst
    flux_err_y = mag_err_y_lsst * flux_y_lsst


    # 2) parameter file

    params = parseParamFile(configfilename, verbose=False, catFilesNeeded=False)

    numB = len(params['bandNames'])
    numObjects = len(galid)

    msg = "get {} objects ".format(numObjects)
    logger.debug(msg)

    logger.debug(params['bandNames'])



    # Generate training data
    #-------------------------


    # what is fluxes and fluxes variance
    fluxes, fluxesVar = np.zeros((numObjects, numB)), np.zeros((numObjects, numB))

    # loop on objects to simulate for the training and save in output training file
    for k in range(numObjects):
        #loop on number of bands
        for i in range(numB):
            if i==0:
                trueFlux = flux_u_lsst[k]
                noise = flux_err_u[k]
            elif i==1:
                trueFlux = flux_g_lsst[k]
                noise = flux_err_g[k]
            elif i==2:
                trueFlux = flux_r_lsst[k]
                noise = flux_err_r[k]
            elif i==3:
                trueFlux = flux_i_lsst[k]
                noise = flux_err_i[k]
            elif i == 4:
                trueFlux = flux_z_lsst[k]
                noise = flux_err_z[k]
            elif i == 5:
                trueFlux = flux_y_lsst[k]
                noise = flux_err_y[k]

            fluxes[k, i] = trueFlux + noise * np.random.randn() # noisy flux
            fluxesVar[k, i] = noise**2.

    # container for training galaxies output
    # at some redshift, provides the flux and its variance inside each band
    data = np.zeros((numObjects, 1 + len(params['training_bandOrder'])))
    bandIndices, bandNames, bandColumns, bandVarColumns, redshiftColumn,refBandColumn = readColumnPositions(params, prefix="training_")

    for ib, pf, pfv in zip(bandIndices, bandColumns, bandVarColumns):
        data[:, pf] = fluxes[:, ib]
        data[:, pfv] = fluxesVar[:, ib]
    data[:, redshiftColumn] = redshifts
    data[:, -1] = 0  # NO type



    msg="write training file {}".format(params['trainingFile'])
    logger.debug(msg)
    np.savetxt(params['trainingFile'], data)





    # Generate Target data : procedure similar to the training
    #-----------------------------------------------------------

    # 1) DESC catalog file
    msg = "read DESC hdf5 file {} ".format(desctargetcatalogfile)
    logger.debug(msg)

    with h5py.File(desctargetcatalogfile, "r") as f:
        # List all groups
        msg = "Keys: %s" % f.keys()
        logger.debug(msg)
        a_group_key = list(f.keys())[0]

        # Get the data
        descdata = list(f[a_group_key])

    logger.debug(descdata)

    f = h5py.File(desctargetcatalogfile, 'r')

    galid = f['photometry/id'][()]
    redshifts = f['photometry/redshift'][()]

    mag_err_g_lsst = f['photometry/mag_err_g_lsst'][()]
    mag_err_i_lsst = f['photometry/mag_err_i_lsst'][()]
    mag_err_r_lsst = f['photometry/mag_err_r_lsst'][()]
    mag_err_u_lsst = f['photometry/mag_err_u_lsst'][()]
    mag_err_y_lsst = f['photometry/mag_err_y_lsst'][()]
    mag_err_z_lsst = f['photometry/mag_err_z_lsst'][()]
    mag_g_lsst = f['photometry/mag_g_lsst'][()]
    mag_i_lsst = f['photometry/mag_i_lsst'][()]
    mag_r_lsst = f['photometry/mag_r_lsst'][()]
    mag_u_lsst = f['photometry/mag_u_lsst'][()]
    mag_y_lsst = f['photometry/mag_y_lsst'][()]
    mag_z_lsst = f['photometry/mag_z_lsst'][()]

    # conversion of magnitude to flux
    flux_g_lsst = np.power(10, -0.4 * mag_g_lsst)
    flux_i_lsst = np.power(10, -0.4 * mag_i_lsst)
    flux_r_lsst = np.power(10, -0.4 * mag_r_lsst)
    flux_u_lsst = np.power(10, -0.4 * mag_u_lsst)
    flux_y_lsst = np.power(10, -0.4 * mag_y_lsst)
    flux_z_lsst = np.power(10, -0.4 * mag_z_lsst)

    # Flux error
    flux_err_u = mag_err_u_lsst * flux_u_lsst
    flux_err_g = mag_err_g_lsst * flux_g_lsst
    flux_err_r = mag_err_r_lsst * flux_r_lsst
    flux_err_i = mag_err_i_lsst * flux_i_lsst
    flux_err_z = mag_err_z_lsst * flux_z_lsst
    flux_err_y = mag_err_y_lsst * flux_y_lsst

    numObjects = len(galid)
    msg = "get {} objects ".format(numObjects)
    logger.debug(msg)

    fluxes, fluxesVar = np.zeros((numObjects, numB)), np.zeros((numObjects, numB))

    # loop on objects in target files
    for k in range(numObjects):
        # loop on bands
        for i in range(numB):
            # compute the flux in that band at the redshift
            if i==0:
                trueFlux = flux_u_lsst[k]
                noise = flux_err_u[k]
            elif i==1:
                trueFlux = flux_g_lsst[k]
                noise = flux_err_g[k]
            elif i==2:
                trueFlux = flux_r_lsst[k]
                noise = flux_err_r[k]
            elif i==3:
                trueFlux = flux_i_lsst[k]
                noise = flux_err_i[k]
            elif i == 4:
                trueFlux = flux_z_lsst[k]
                noise = flux_err_z[k]
            elif i == 5:
                trueFlux = flux_y_lsst[k]
                noise = flux_err_y[k]


            fluxes[k, i] = trueFlux + noise * np.random.randn()
            fluxesVar[k, i] = noise**2.

    data = np.zeros((numObjects, 1 + len(params['target_bandOrder'])))
    bandIndices, bandNames, bandColumns, bandVarColumns, redshiftColumn,refBandColumn = readColumnPositions(params, prefix="target_")

    for ib, pf, pfv in zip(bandIndices, bandColumns, bandVarColumns):
        data[:, pf] = fluxes[:, ib]
        data[:, pfv] = fluxesVar[:, ib]
    data[:, redshiftColumn] = redshifts
    data[:, -1] = 0  # NO TYPE

    msg = "write file {}".format(os.path.basename(params['targetFile']))
    logger.debug(msg)

    msg = "write target file {}".format(params['targetFile'])
    logger.debug(msg)

    np.savetxt(params['targetFile'], data)


if __name__ == "__main__":
    # execute only if run as a script


    msg="Start convertDESCcat.py"
    logger.info(msg)
    logger.info("--- convert DESC catalogs ---")



    if len(sys.argv) < 4:
        raise Exception('Please provide a parameter file and the tranning and catalog files')

    convertDESCcat(sys.argv[1],sys.argv[2],sys.argv[3])
