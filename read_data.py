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

from astropy.table import Table

import fitting

co_filename = __file__
co_path     = os.path.dirname(co_filename) + '/'

def OH_flag_arr(OH_file, x0, silent=False, verbose=False):
    '''
    Function to get information on OH night skylines

    Parameters
    ----------
    OH_file : string
      Filename containing wavelengths affected by OH night skylines

    x0 : numpy array
      Wavelength grid in Angstrom

    Returns
    -------
    OH_dict0 : dict
      Dictionary of arrays containing OH skyline info:
        'OH_flag0' - Array containing 1 and 0. Same dimensions as [x0]
        'OH_xmin0' - Array containing minimum wavelength for night skyline
        'OH_xmax0' - Array containing maximum wavelength for night skyline

    Notes
    -----
    Created by Chun Ly, 10 February 2017
    '''

    if silent == True: print '## Reading : ', OH_file
    OH_data = asc.read(OH_file)

    OH_xmin0 = OH_data['col1'].data.tolist()
    OH_xmax0 = OH_data['col2'].data.tolist()

    if silent == False:
        print '## Adding in A- and B-band of atmospheric absorption'
        OH_xmin0 += [7590.0, 6860.0]
        OH_xmax0 += [7670.0, 6890.0]

    OH_flag0 = np.zeros(len(x0), dtype=np.int8)
    for oo in range(len(OH_xmin0)):
        mark2 = np.where((x0 >= OH_xmin0[oo]) & (x0 <= OH_xmax0[oo]))[0]
        if len(mark2) > 0: OH_flag[mark2] = 1

    OH_dict0 = {'OH_flag0':OH_flag0, 'OH_xmin0':OH_xmin0, 'OH_xmax0':OH_xmax0}
    return OH_dict0
#enddef

def get_tagnames(resol='low', DEEP2=False, SDF=True, weak=False,
                 silent=False, verbose=True):
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
    Created by Chun Ly, 9 February 2017
    '''

    # Contains information about fitted emission lines
    # Later + on 09/02/2017
    line_fit_file = co_path + 'fit_lines.'+resol+'res.txt'
    if silent == False: print '### Reading : ', line_fit_file
    fit_data0 = asc.read(line_fit_file, format='commented_header')

    if verbose == True: print fit_data0

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

    return tagnames, dtype, fit_data0
#enddef

def main(infile0, cat_file0, zspec0=None, OH_file=None, resol='low',
         out_pdf=None, DEEP2=False, SDF=True, weak=False, silent=False,
         verbose=True):

    '''
    Main function for reading in data

    Parameters
    ----------
    infile0 : string
      Filename for FITS 2-D image that contains 1-D spectra

    cat_file0 : string
      Filename for FITS binary table with individual entries corresponding
      to each 1-D spectra in infile0

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
    Modified by Chun Ly, 9 February 2017
     - Add cat_file0 input
     - Fill in [emline_data] for basic identification columns
     - Added out_pdf keyword option
    Modified by Chun Ly, 10 February 2017
     - Add spec_file to dict0 for fitting.main()
     - Add call to OH_flag_arr()
    '''

    if silent == False: print '### Begin read_data.main | '+systime()

    if silent == False: print '### Reading : ', infile0
    data0, hdr0 = fits.getdata(infile0, header=True)

    # + on 09/02/2017
    if silent == False: print '### Reading : ', cat_file0
    cat_data0, cat_hdr0 = fits.getdata(cat_file0, header=True)
    cat_data0 = Table(cat_data0)

    if zspec0 == None: zspec0 = np.zeros(len(data0))

    # + on 09/02/2017
    if out_pdf == None:
        str_in = '.fits.gz' if '.gz' in infile0 else '.fits'
        out_pdf = infile0.replace(str_in,'.pdf')
    if silent == False: print '### out_pdf : ', out_pdf

    l0    = hdr0['CRVAL1'] # Wavelength minimum
    dl    = hdr0['CDELT1'] # Wavelength dispersion
    n_pix = hdr0['NAXIS1'] # Number of pixels
    if silent == False:
        print '## Minimum wavelength : ', l0
        print '## Wavelength dispersion : ', dl
        print '## Number of spectral pixels : ', n_pix

    x0 = l0 + dl * np.arange(n_pix)

    # Get OH skyline information | + on 10/02/2017
    if OH_file != None: OH_dict0 = OH_flag_arr(OH_file, x0)

    n_spec = len(data0)

    tagnames, dtype, \
        fit_data0 = get_tagnames(resol=resol, DEEP2=DEEP2, SDF=SDF,
                                 weak=weak, silent=silent, verbose=verbose)

    # Define columns of data and table
    arr0 = {}
    for tt in range(len(tagnames)):
        if dtype[tt] == 'i':
            arr0[tagnames[tt]] = np.zeros(n_spec, dtype=np.int)
        if 'S' in dtype[tt]:
            arr0[tagnames[tt]] = np.repeat('N/A       ',n_spec)
        if 'f' in dtype[tt] or 'e' in dtype[tt]:
            arr0[tagnames[tt]] = np.zeros(n_spec, dtype=np.float32)

    emline_data = Table(arr0, names=tagnames)

    # + on 09/02/2017
    init_tags = ['OBJNO', 'AP', 'SLIT']
    for tag in init_tags:
        if tag in tagnames: emline_data[tag] = cat_data0[tag]

    emline_data['ZSPEC'] = zspec0

    dict0 = {'data0':data0, 'emline_data':emline_data, 'fit_data0':fit_data0,
             'spec_file': os.path.basename(infile0), 'x0': x0}

    # + on 10/02/2017
    if OH_file != None: dict0['OH_dict0'] = OH_dict0

    emline_data = fitting.main(dict0, out_pdf, silent=silent, verbose=verbose)

    if silent == False: print '### End read_data.main | '+systime()
#enddef

