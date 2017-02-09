"""
read_data
=========

Set of functions to read in FITS 1-D data
"""

import sys, os

from chun_codes import systime

from os.path import exists
import commands
from astropy.io import ascii as asc
from astropy.io import fits

import numpy as np
import array, time, sets

import matplotlib.pyplot as plt
import glob

co_filename = __file__
co_path     = os.path.dirname(co_filename) + '/'

def get_tagnames(resol='low', DEEP2=False, SDF=True, weak=False):
    '''
    Function to get tagnames for table with emission-line fitting results

    Parameters
    ----------
    resol : string
      Spectral resolution of 1-D spectra. For low resolution, the
      [OII] doublet cannot be resolved. Default: 'low'

    DEEP2 : bool
      Indicate whether data is for DEEP2. Default: False

    SDF : bool
      Indicate whether data is for the SDF. Default: True

    weak : bool
      Indicate whether to fit weak nebular emission lines. Default: False

    Returns
    -------
    tagnames : list
      List of table column names

    dtype : list
      List of data types for table column names

    Notes
    -----
    Created by Chun Ly, 8 February 2017
    '''

    if DEEP2 == True:
        tagnames = ['OBJNO', 'SLIT', 'LINE', 'ZSPEC']
        dtype    = ['i', 'S10', 'i', 'f8']
    if SDF == True:
        tagnames = ['AP', 'LINE', 'ZSPEC']
        dtype    = ['S10', 'i', 'f8']

    cols    = ['LAMBDA', 'ZSPEC', 'PEAK', 'SIGMA', 'Y0',
               'FLUX', 'FLUX_ERR', 'FLUX_DATA', 'NOISE', 'SNR']
    c_dtype = ['f8','f8','e','f8','e','e','e','e','e','f8']

    if resol == 'high':
        t_tags   = ['OII_3727_'+str0 for str0 in cols]
        tagnames = tagnames + t_tags
        dtype    = dtype + c_dtype

    for ll in range(len(fit_data0)):
        t_tags   = [fit_data0['lines_txt'][ll]+'_'+str0 for str0 in cols]
        tagnames = tagnames + t_tags
        dtype    = dtype + c_dtype

    if weak == False:
        b_tags = ['H8_ABS_EW', 'H7_ABS_EW', 'H_DELTA_ABS_EW',
                  'H_GAMMA_ABS_EW', 'H_BETA_ABS_EW', 'H_ALPHA_ABS_EW']
        b_dtype = ['f8', 'f8', 'f8', 'f8', 'f8', 'f8']
    else:
        b_tags = ['H10_ABS_EW', 'H9_ABS_EW']
        b_dtype = ['f8', 'f8']

    tagnames = tagnames + b_tags
    dtype    = dtype + b_dtype

    # Add columns for spectral coverage min/max wavelength
    wave_tags = ['LMIN', 'LMAX', 'LMIN0', 'LMAX0']
    w_dtype   = ['f8', 'f8', 'f8', 'f8']

    tagnames  = tagnames + wave_tags
    dtype     = dtype + w_dtype

    return tagnames, dtype
#enddef

def main(infile0, zspec0=None, OH_file=None, resol='low', DEEP2=False,
         SDF=True, weak=False, silent=False, verbose=True):

    '''
    Main function for reading in data

    Parameters
    ----------
    infile0 : string
      Filename for FITS data that contains 1-D spectra

    zspec0 : list or numpy array
      Array of spectroscopic redshifts for each 1-D spectrum
      Default: np.array of zeros

    OH_file : string
      Filename containing wavelengths affected by OH night skylines

    resol : string
      Spectral resolution of 1-D spectra. For low resolution, the
      [OII] doublet cannot be resolved. Default: 'low'

    DEEP2 : bool
      Indicate whether data is for DEEP2. Default: False

    SDF : bool
      Indicate whether data is for the SDF. Default: True

    weak : bool
      Indicate whether to fit weak nebular emission lines. Default: False

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 8 February 2017
    '''

    if silent == False: print '### Begin read_data.main | '+systime()

    if silent == False: print '### Reading : ', infile0
    data0, hdr0 = fits.getdata(infile0, header=True)

    if zspec0 == None: zspec0 = np.zeros(len(data0))

    l0    = hdr0['CRVAL1'] # Wavelength minimum
    dl    = hdr0['CDELT1'] # Wavelength dispersion
    n_pix = hdr0['NAXIS1'] # Number of pixels
    if silent == False:
        print '## Minimum wavelength : ', l0
        print '## Wavelength dispersion : ', dl
        print '## Number of spectral pixels : ', n_pix

    x0 = l0 + dl * np.arange(n_pix)

    n_spec = len(data0)

    line_fit_file = co_path + 'fit_lines.'+resol+'res.txt'
    if silent == False: print '### Reading : ', line_fit_file
    fit_data0 = asc.read(line_fit_file, format='commented_header')

    if verbose == True: print fit_data0

    tagnames, dtype = get_tagnames(resol=resol, DEEP2=DEEP2, SDF=SDF,
                                   weak=weak)

    # Define columns of data and table
    arr0 = {}
    for tt in range(len(tagnames)):
        if dtype[tt] == 'i':
            arr0[tagnames[tt]] = np.zeros(n_spec, dtype=np.int)
        if 'S' in dtype[tt]:
            arr0[tagnames[tt]] = np.repeat('N/A       ',n_spec)
        if 'f' in dtype[tt] or 'e' in dtype[tt]:
            arr0[tagnames[tt]] = np.zeros(n_spec, dtype=np.float32)

    emline_data = Table(arr0)
    
    if silent == False: print '### End read_data.main | '+systime()
#enddef

