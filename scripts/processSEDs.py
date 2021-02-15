####################################################################################################
#
# script : processSED.py
#
# process the library of SEDs and project them onto the filters, (for the mean fct of the GP)
# -  for each readhift in the grid and for for each filter selected in config
# (which may take a few minutes depending on the settings you set):
#
# output file : sed_name + '_fluxredshiftmod.txt'
######################################################################################################

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


from delight.io import *
from delight.utils import *

logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG', logger=logger,fmt='%(asctime)s,%(msecs)03d %(programname)s, %(name)s[%(process)d] %(levelname)s %(message)s')


if len(sys.argv) < 2:
    raise Exception('Please provide a parameter file')

logger.info("--- Process SED ---")

params = parseParamFile(sys.argv[1], verbose=False, catFilesNeeded=False)
bandNames = params['bandNames']
dir_seds = params['templates_directory']
dir_filters = params['bands_directory']
lambdaRef = params['lambdaRef']
sed_names = params['templates_names']
fmt = '.dat'

DL = approx_DL()

#redshift grid
redshiftDistGrid, redshiftGrid, redshiftGridGP = createGrids(params)
numZ = redshiftGrid.size

# Loop over SEDs
# create a file per SED of all possible flux in band
for sed_name in sed_names:
    seddata = np.genfromtxt(dir_seds + '/' + sed_name + fmt)
    seddata[:, 1] *= seddata[:, 0] # SDC : multiply luminosity by wl ?
    # SDC: OK if luminosity is in wl bins ! To be checked !!!!
    ref = np.interp(lambdaRef, seddata[:, 0], seddata[:, 1])
    seddata[:, 1] /= ref  # normalisation at lambdaRef
    sed_interp = interp1d(seddata[:, 0], seddata[:, 1]) # interpolation function

    # container for the model flux at redshifts in grid
    # each row correspond to fluxes in the different bands at a  fixed redshift
    # redshift along row, fluxes along column
    f_mod = np.zeros((redshiftGrid.size, len(bandNames)))

    # Loop over bands
    # jf index on bands
    for jf, band in enumerate(bandNames):
        fname_in = dir_filters + '/' + band + '.res'
        data = np.genfromtxt(fname_in)
        xf, yf = data[:, 0], data[:, 1]  # filter transmission function at observation
        #yf /= xf  # divide by lambda
        # Only consider range where >1% max
        ind = np.where(yf > 0.01*np.max(yf))[0]
        lambdaMin, lambdaMax = xf[ind[0]], xf[ind[-1]]
        norm = np.trapz(yf/xf, x=xf) # SDC: probably Cb by integrating with trapeze method

        ###############################
        # iz index on redshift
        # hope no lambda factor is forgetted
        ########################
        for iz in range(redshiftGrid.size):
            opz = (redshiftGrid[iz] + 1)  #(1+z) factor
            # wavelength in lambda_em
            xf_z = np.linspace(lambdaMin / opz, lambdaMax / opz, num=5000)
            yf_z = interp1d(xf / opz, yf)(xf_z) # filter interpol function in lambda_em args
            ysed = sed_interp(xf_z) # sed
            f_mod[iz, jf] = np.trapz(ysed * yf_z, x=xf_z) / norm  # trapez integration of lum over filter
            f_mod[iz, jf] *= opz**2. / DL(redshiftGrid[iz])**2. / (4*np.pi) # output if flux units (AB system ?)
    # for each SED, save the flux at each redshift (along row) for each
    np.savetxt(dir_seds + '/' + sed_name + '_fluxredshiftmod.txt', f_mod)
