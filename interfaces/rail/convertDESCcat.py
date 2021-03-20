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

from interfaces.rail.utils  import load_training_data, get_input_data_size_hdf5,load_raw_hdf5_data

logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG', logger=logger,fmt='%(asctime)s,%(msecs)03d %(programname)s, %(name)s[%(process)d] %(levelname)s %(message)s')


def group_entries(f):
    """
    group entries in single numpy array

    """
    galid = f['id'][()][:, np.newaxis]
    redshift = f['redshift'][()][:, np.newaxis]
    mag_err_g_lsst = f['mag_err_g_lsst'][()][:, np.newaxis]
    mag_err_i_lsst = f['mag_err_i_lsst'][()][:, np.newaxis]
    mag_err_r_lsst = f['mag_err_r_lsst'][()][:, np.newaxis]
    mag_err_u_lsst = f['mag_err_u_lsst'][()][:, np.newaxis]
    mag_err_y_lsst = f['mag_err_y_lsst'][()][:, np.newaxis]
    mag_err_z_lsst = f['mag_err_z_lsst'][()][:, np.newaxis]
    mag_g_lsst = f['mag_g_lsst'][()][:, np.newaxis]
    mag_i_lsst = f['mag_i_lsst'][()][:, np.newaxis]
    mag_r_lsst = f['mag_r_lsst'][()][:, np.newaxis]
    mag_u_lsst = f['mag_u_lsst'][()][:, np.newaxis]
    mag_y_lsst = f['mag_y_lsst'][()][:, np.newaxis]
    mag_z_lsst = f['mag_z_lsst'][()][:, np.newaxis]

    full_arr = np.hstack((galid, redshift, mag_u_lsst, mag_g_lsst, mag_r_lsst, mag_i_lsst, mag_z_lsst, mag_y_lsst, \
                          mag_err_u_lsst, mag_err_g_lsst, mag_err_r_lsst, mag_err_i_lsst, mag_err_z_lsst,
                          mag_err_y_lsst))
    return full_arr


def filter_mag_entries(d,nb=6):
    """
    Filter bad data with bad magnitudes

    input
      - d: array of magnitudes and errors
      - nb : number of bands
    output :
      - indexes of row to be filtered

    """

    u = d[:, 2]
    idx_u = np.where(u > 31.8)[0]

    return idx_u


def mag_to_flux(d,nb=6):
    """

    Convert magnitudes to fluxes

    input:
       -d : array of magnitudes with errors


    :return:
        array of fluxes with error
    """

    fluxes = np.zeros_like(d)

    fluxes[:, 0] = d[:, 0]  # object index
    fluxes[:, 1] = d[:, 1]  # redshift

    for idx in np.arange(nb):
        fluxes[:, 2 + idx] = np.power(10, -0.4 * d[:, 2 + idx]) # fluxes
        fluxes[:, 8 + idx] = fluxes[:, 2 + idx] * d[:, 8 + idx] # errors on fluxes
    return fluxes



def filter_flux_entries(d,nb=6,nsig=5):
    """
    Filter noisy data on the the number SNR

    input :
     - d: flux and errors array
     - nb : number of bands
     - nsig : number of sigma

     output:
        indexes of row to suppress

    """


    # collection of indexes
    indexes = []
    indexes = np.array(indexes, dtype=np.int)

    for idx in np.arange(nb):
        ratio = d[:, 2 + idx] / d[:, 8 + idx]  # flux divided by sigma-flux
        bad_indexes = np.where(ratio < nsig)[0]
        indexes = np.concatenate((indexes, bad_indexes))

    indexes = np.unique(indexes)
    return np.sort(indexes)


def convertDESCcatChunk(configfilename,data,chunknum):

        """

        Convert files in ascii format to be used by Delight

        input args:
        - configfilename : Delight configuration file containg path for output files (flux variances and redshifts)
        - data
        - chunknum : number of the chunk

        output :
        the Delight training and target file which path is in configuration file

        :param configfilename:
        :return:
        """
        msg="--- Convert DESC catalogs chunck {}---".format(chunknum)
        logger.info(msg)


        # produce a numpy array
        magdata = group_entries(data)
        # filter bad data
        indexes_bad = filter_mag_entries(magdata)
        magdata_f = np.delete(magdata, indexes_bad, axis=0)

        # convert mag to fluxes
        fdata = mag_to_flux(magdata_f)

        # filter bad data
        indexes_bad = filter_flux_entries(fdata)

        magdata_f = np.delete(magdata_f, indexes_bad, axis=0)
        fdata_f = np.delete(fdata, indexes_bad, axis=0)

        gid = magdata_f[:, 0]
        rs = magdata_f[:, 1]

        # 2) parameter file

        params = parseParamFile(configfilename, verbose=False, catFilesNeeded=False)

        numB = len(params['bandNames'])
        numObjects = len(gid)

        msg = "get {} objects ".format(numObjects)
        logger.debug(msg)

        logger.debug(params['bandNames'])

        # Generate target data
        # -------------------------

        # what is fluxes and fluxes variance
        fluxes, fluxesVar = np.zeros((numObjects, numB)), np.zeros((numObjects, numB))

        # loop on objects to simulate for the target and save in output trarget file
        for k in range(numObjects):
            # loop on number of bands
            for i in range(numB):
                trueFlux = fdata_f[k, 2 + i]
                noise = fdata_f[k, 8 + i]

                # fluxes[k, i] = trueFlux + noise * np.random.randn() # noisy flux
                fluxes[k, i] = trueFlux

                if fluxes[k, i] < 0:
                    # fluxes[k, i]=np.abs(noise)/10.
                    fluxes[k, i] = trueFlux

                fluxesVar[k, i] = noise ** 2.

        # container for target galaxies output
        # at some redshift, provides the flux and its variance inside each band
        

        data = np.zeros((numObjects, 1 + len(params['target_bandOrder'])))
        bandIndices, bandNames, bandColumns, bandVarColumns, redshiftColumn, refBandColumn = readColumnPositions(params,
                                                                                                                 prefix="target_")

        for ib, pf, pfv in zip(bandIndices, bandColumns, bandVarColumns):
            data[:, pf] = fluxes[:, ib]
            data[:, pfv] = fluxesVar[:, ib]
        data[:, redshiftColumn] = rs
        data[:, -1] = 0  # NO TYPE

        msg = "write file {}".format(os.path.basename(params['targetFile']))
        logger.debug(msg)

        msg = "write target file {}".format(params['targetFile'])
        logger.debug(msg)

        outputdir = os.path.dirname(params['targetFile'])
        if not os.path.exists(outputdir):
            msg = " outputdir not existing {} then create it ".format(outputdir)
            logger.info(msg)
            os.makedirs(outputdir)

        np.savetxt(params['targetFile'], data)



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

    f = load_raw_hdf5_data(desctraincatalogfile, groupname='photometry')

    # produce a numpy array
    magdata = group_entries(f)
    # filter bad data
    indexes_bad = filter_mag_entries(magdata)
    magdata_f = np.delete(magdata, indexes_bad, axis=0)

    # convert mag to fluxes
    fdata = mag_to_flux(magdata_f)

    #filter bad data
    indexes_bad = filter_flux_entries(fdata)

    magdata_f = np.delete(magdata_f, indexes_bad, axis=0)
    fdata_f = np.delete(fdata, indexes_bad, axis=0)


    gid = magdata_f[:, 0]
    rs = magdata_f[:, 1]




    # 2) parameter file

    params = parseParamFile(configfilename, verbose=False, catFilesNeeded=False)

    numB = len(params['bandNames'])
    numObjects = len(gid)

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
            trueFlux = fdata_f[k,2+i]
            noise    = fdata_f[k,8+i]

            #fluxes[k, i] = trueFlux + noise * np.random.randn() # noisy flux
            fluxes[k, i] = trueFlux

            if fluxes[k, i]<0:
                #fluxes[k, i]=np.abs(noise)/10.
                fluxes[k, i] = trueFlux

            fluxesVar[k, i] = noise**2.

    # container for training galaxies output
    # at some redshift, provides the flux and its variance inside each band
    data = np.zeros((numObjects, 1 + len(params['training_bandOrder'])))
    bandIndices, bandNames, bandColumns, bandVarColumns, redshiftColumn,refBandColumn = readColumnPositions(params, prefix="training_")

    for ib, pf, pfv in zip(bandIndices, bandColumns, bandVarColumns):
        data[:, pf] = fluxes[:, ib]
        data[:, pfv] = fluxesVar[:, ib]
    data[:, redshiftColumn] = rs
    data[:, -1] = 0  # NO type


    msg="write training file {}".format(params['trainingFile'])
    logger.debug(msg)

    outputdir=os.path.dirname(params['trainingFile'])
    if not os.path.exists(outputdir):
        msg = " outputdir not existing {} then create it ".format(outputdir)
        logger.info(msg)
        os.makedirs(outputdir)


    np.savetxt(params['trainingFile'], data)




    # Generate Target data : procedure similar to the training
    #-----------------------------------------------------------

    # 1) DESC catalog file
    msg = "read DESC hdf5 file {} ".format(desctargetcatalogfile)
    logger.debug(msg)

    f = load_raw_hdf5_data(desctargetcatalogfile, groupname='photometry')

    # produce a numpy array
    magdata = group_entries(f)
    # filter bad data
    indexes_bad = filter_mag_entries(magdata)
    magdata_f = np.delete(magdata, indexes_bad, axis=0)

    # convert mag to fluxes
    fdata = mag_to_flux(magdata_f)

    # filter bad data
    indexes_bad = filter_flux_entries(fdata)

    magdata_f = np.delete(magdata_f, indexes_bad, axis=0)
    fdata_f = np.delete(fdata, indexes_bad, axis=0)

    gid = magdata_f[:, 0]
    rs = magdata_f[:, 1]





    numObjects = len(gid)
    msg = "get {} objects ".format(numObjects)
    logger.debug(msg)

    fluxes, fluxesVar = np.zeros((numObjects, numB)), np.zeros((numObjects, numB))

    # loop on objects in target files
    for k in range(numObjects):
        # loop on bands
        for i in range(numB):
            # compute the flux in that band at the redshift
            trueFlux = fdata_f[k, 2 + i]
            noise = fdata_f[k, 8 + i]

            #fluxes[k, i] = trueFlux + noise * np.random.randn()
            fluxes[k, i] = trueFlux

            if fluxes[k, i]<0:
                #fluxes[k, i]=np.abs(noise)/10.
                fluxes[k, i] = trueFlux

            fluxesVar[k, i] = noise**2.

    data = np.zeros((numObjects, 1 + len(params['target_bandOrder'])))
    bandIndices, bandNames, bandColumns, bandVarColumns, redshiftColumn,refBandColumn = readColumnPositions(params, prefix="target_")

    for ib, pf, pfv in zip(bandIndices, bandColumns, bandVarColumns):
        data[:, pf] = fluxes[:, ib]
        data[:, pfv] = fluxesVar[:, ib]
    data[:, redshiftColumn] = rs
    data[:, -1] = 0  # NO TYPE

    msg = "write file {}".format(os.path.basename(params['targetFile']))
    logger.debug(msg)

    msg = "write target file {}".format(params['targetFile'])
    logger.debug(msg)

    outputdir = os.path.dirname(params['targetFile'])
    if not os.path.exists(outputdir):
        msg = " outputdir not existing {} then create it ".format(outputdir)
        logger.info(msg)
        os.makedirs(outputdir)

    np.savetxt(params['targetFile'], data)


if __name__ == "__main__":
    # execute only if run as a script


    msg="Start convertDESCcat.py"
    logger.info(msg)
    logger.info("--- convert DESC catalogs ---")



    if len(sys.argv) < 4:
        raise Exception('Please provide a parameter file and the training and validation and catalog files')

    convertDESCcat(sys.argv[1],sys.argv[2],sys.argv[3])
