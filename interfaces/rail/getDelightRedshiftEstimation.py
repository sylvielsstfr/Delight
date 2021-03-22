import sys
import os
import numpy as np
from functools import reduce

import pprint

from delight.io import *
from delight.utils import *
import h5py


logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG', logger=logger,fmt='%(asctime)s,%(msecs)03d %(programname)s, %(name)s[%(process)d] %(levelname)s %(message)s')



def getDelightRedshiftEstimation(configfilename,chunknum,nsize,index_sel):
    """
    zmode,widths = getDelightRedshiftEstimation(delightparamfilechunk,self.chunknum,d,indexes_sel)

    :return:
    """

    msg = "--- getDelightRedshiftEstimation({}) for chunk {}---".format(nsize,chunknum)
    logger.info(msg)

    zmode = np.full(nsize, fill_value=-1,dtype=np.float)

    params = parseParamFile(configfilename, verbose=False)

    redshiftDistGrid, redshiftGrid, redshiftGridGP = createGrids(params)

    pdfs = np.loadtxt(params['redshiftpdfFile'])
    pdfs /= np.trapz(pdfs, x=redshiftGrid, axis=1)[:, None]


    loc=-1
    for idx in np.arange(nsize):
        if idx in index_sel:
            loc+=1
            idx_zmode = np.where( pdfs[loc] == pdfs[loc].max() )[0][0]

            if idx_zmode<len(redshiftGrid):
                thez=redshiftGrid[idx_zmode]
            else:
                thez = -1

            zmode[idx] = thez



    widths = 0.1 * (1.0 + zmode)
    widths =np.where(zmode>0.0,widths,-1)

    return zmode,widths

