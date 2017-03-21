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
from astropy import log # + on 01/03/2017

import numpy as np

import matplotlib.pyplot as plt
import glob

from astropy.table import Table

import fitting

from . import cols

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

    if silent == True: log.info('## Reading : '+OH_file)
    OH_data = asc.read(OH_file)

    OH_xmin0 = OH_data['col1'].data.tolist()
    OH_xmax0 = OH_data['col2'].data.tolist()

    if silent == False:
        log.info('## Adding in A- and B-band of atmospheric absorption')
    OH_xmin0 += [7590.0, 6860.0]
    OH_xmax0 += [7670.0, 6890.0]

    OH_flag0 = np.zeros(len(x0), dtype=np.int8)
    for oo in range(len(OH_xmin0)):
        mark2 = np.where((x0 >= OH_xmin0[oo]) & (x0 <= OH_xmax0[oo]))[0]
        if len(mark2) > 0: OH_flag0[mark2] = 1

    OH_dict0 = {'OH_flag0':OH_flag0, 'OH_xmin0':OH_xmin0, 'OH_xmax0':OH_xmax0}
    return OH_dict0
#enddef

def get_tagnames(init_dict0, resol='low', weak=False, silent=False,
                 verbose=True):
    '''
    Function to get tagnames for table with emission-line fitting results

    Parameters
    ----------
    resol : string
      Spectral resolution of 1-D spectra. For low resolution, the
      [OII] doublet cannot be resolved. Default: 'low'

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
    Modified by Chun Ly, 11 February 2017
     - Add init_dict0 to allow for user-defined arrays to pass forward
    Modified by Chun Ly, 12 February 2017
     - Added line_type to indicate which columns are associated with each
       emission line
    '''

    # Contains information about fitted emission lines
    # Later + on 09/02/2017
    line_fit_file = co_path + 'fit_lines.'+resol+'res.txt'
    if silent == False: log.info('### Reading : '+line_fit_file)
    fit_data0 = asc.read(line_fit_file, format='commented_header')

    if verbose == True: print fit_data0

    # Mod on 11/02/2107. Switched to dtype as first dict key
    tagnames  = init_dict0.keys()[1:]
    dtype     = init_dict0['dtype']
    line_type = ['N/A']*len(tagnames)
    #if DEEP2 == True:
    #    tagnames = ['OBJNO', 'SLIT', 'LINE', 'ZSPEC']
    #    dtype    = ['i', 'S10', 'i', 'f8']
    #if SDF == True:
    #    tagnames = ['AP', 'LINE', 'ZSPEC']
    #    dtype    = ['S10', 'i', 'f8']

    # cols    = ['LAMBDA', 'ZSPEC', 'PEAK', 'SIGMA', 'Y0',
    #           'FLUX', 'FLUX_ERR', 'FLUX_DATA', 'NOISE', 'SNR']
    c_dtype = ['f8','f8','e','f8','e','e','e','e','e','f8']

    if resol == 'high':
        t_tags     = ['OII_3727_'+str0 for str0 in cols]
        tagnames   = tagnames + t_tags
        dtype      = dtype + c_dtype
        line_type += ['OII_3727'] * len(t_tags)

    for ll in range(len(fit_data0)):
        t_tags     = [fit_data0['lines_txt'][ll]+'_'+str0 for str0 in cols]
        tagnames   = tagnames + t_tags
        dtype      = dtype + c_dtype
        line_type += [fit_data0['lines_txt'][ll]] * len(t_tags)

    if weak == False:
        b_tags = ['H8_ABS_EW', 'H7_ABS_EW', 'H_DELTA_ABS_EW',
                  'H_GAMMA_ABS_EW', 'H_BETA_ABS_EW', 'H_ALPHA_ABS_EW']
        b_dtype = ['f8', 'f8', 'f8', 'f8', 'f8', 'f8']
    else:
        b_tags = ['H10_ABS_EW', 'H9_ABS_EW']
        b_dtype = ['f8', 'f8']

    tagnames   = tagnames + b_tags
    dtype      = dtype + b_dtype
    line_type += ['N/A'] * len(b_tags)

    # Add columns for spectral coverage min/max wavelength
    wave_tags = ['LMIN', 'LMAX', 'LMIN0', 'LMAX0']
    w_dtype   = ['f8', 'f8', 'f8', 'f8']

    tagnames   = tagnames + wave_tags
    dtype      = dtype + w_dtype
    line_type += ['N/A'] * len(wave_tags)

    return tagnames, dtype, fit_data0, line_type
#enddef

def main(infile0, init_dict0, OH_file=None, resol='low', out_pdf=None,
         out_fits=None, weak=False, silent=False, verbose=True):

    '''
    Main function for reading in data

    Parameters
    ----------
    infile0 : string
      Filename for FITS 2-D image that contains 1-D spectra

    init_dict0 : collections.OrderedDict()
      Dictionary containing arrays to pass into emission-line table.
      This can include a source name/ID, RA, Dec, or other
      identification (e.g., SLIT). It must include 'ZSPEC', and
      'dtype'. dtype indicates what data format for each column in the
      same order to initialize the values for the final emission-line
      table.
        - ZSPEC must be provided at high precision, to four decimal places
          to avoid fitting nearby OH skylines. This code uses the zspec
          as a prior. If not provided, it is assumed that the spectra
          have been shifted to the rest-frame and all zspec values will be
          ZERO.
        - dtype must be placed at the *beginning* of the OrderedDict().
        - dtype must be a list for concatenation purposes
        - 'LINE' is an index array providing which spectra to perform the
          fitting. The indexing should begin from 1 (for the first spectra)
          and continue on. If it is not provided, it will assume that all
          spectra will be fitted.
        - Note that the input arrays must have the same dimensions. However,
          this can be different from the size of the FITS 2-D image (infile0)
      Below is an example:
        init_dict0 = collections.OrderedDict()
        dtype = ['i', 'S8', 'i', 'f8']
        init_dict0['dtype'] = dtype
        init_dict0['ID']    = ID
        init_dict0['SLIT']  = slit
        init_dict0['LINE']  = line
        init_dict0['ZSPEC'] = zspec

    OH_file : string
      Filename containing wavelengths affected by OH night skylines

    resol : string
      Spectral resolution of 1-D spectra. For low resolution, the
      [OII] doublet cannot be resolved. Default: 'low'

    out_pdf : string
      Filename for output PDF file

    out_fits : string
      Filename for output FITS file containing emission-line fitting results

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
    Modified by Chun Ly, 11 February 2017
     - Added init_dict0 input for user-defined arrays to pass to
       emission-line table
     - Delete DEEP2 and SDF keyword options
    Modified by Chun Ly, 12 February 2017
     - Handle case when 'LINE' is not provided in init_dict0
     - Handle case when 'ZSPEC' is not provided in init_dict0 -- Assume rest-frame
     - Add out_fits keyword option
    Modified by Chun Ly, 2 March 2017
     - Include [dl] in dict0
    '''

    if silent == False: log.info('### Begin read_data.main : '+systime())

    # + on 11/02/2017
    if 'dtype' not in init_dict0.keys():
        log.error('## dtype not provided in init_dict0!!!')
        log.error('## Exiting!!!')
        return

    if silent == False: log.info('### Reading : '+infile0)
    data0, hdr0 = fits.getdata(infile0, header=True)

    # + on 12/02/2017
    if 'LINE' not in init_dict0.keys():
        if silent == False: log.info('### Adding default LINE index')
        init_dict0['LINE']  = 1+np.arange(len(data0))
        init_dict0['dtype'] = init_dict0['dtype'] + ['i']
        if silent == False: log.info(init_dict0.keys())
    line = init_dict0['LINE']

    # Mod on 12/02/2017
    if 'ZSPEC' not in init_dict0.keys():
        init_dict0['ZSPEC'] = np.zeros(len(line))

    str_in  = '.fits.gz' if '.gz' in infile0 else '.fits'
    # + on 09/02/2017
    if out_pdf == None: out_pdf = infile0.replace(str_in,'.ELF.pdf')
    if silent == False: log.info('### out_pdf : '+out_pdf)

    # + on 12/02/2017
    if out_fits == None: out_fits = infile0.replace(str_in,'.ELF.fits')
    if silent == False: log.info('### out_fits : '+out_fits)

    l0    = hdr0['CRVAL1'] # Wavelength minimum
    dl    = hdr0['CDELT1'] # Wavelength dispersion
    n_pix = hdr0['NAXIS1'] # Number of pixels
    if silent == False:
        log.info('## Minimum wavelength : %f' % l0)
        log.info('## Wavelength dispersion : %f' % dl)
        log.info('## Number of spectral pixels : %i' % n_pix)

    x0 = l0 + dl * np.arange(n_pix)

    # Get OH skyline information | + on 10/02/2017
    if OH_file != None: OH_dict0 = OH_flag_arr(OH_file, x0)

    n_line = len(line) #len(data0)

    # Mod on 11/02/2017, 12/02/2017
    tagnames, dtype, fit_data0, \
        line_type = get_tagnames(init_dict0, resol=resol, weak=weak,
                                 silent=silent, verbose=verbose)

    # Define columns of data and table
    arr0 = {}
    for tt in range(len(tagnames)):
        if dtype[tt] == 'i':
            arr0[tagnames[tt]] = np.zeros(n_line, dtype=np.int)
        if 'S' in dtype[tt]:
            arr0[tagnames[tt]] = np.repeat('N/A       ',n_line)
        if 'f' in dtype[tt] or 'e' in dtype[tt]:
            arr0[tagnames[tt]] = np.zeros(n_line, dtype=np.float32)

    emline_data = Table(arr0, names=tagnames)

    # + on 09/02/2017. Mod on 11/02/2017
    init_tags = init_dict0.keys()[1:]
    if silent == False: log.info('### Filling in : '+", ".join(init_tags))
    for tag in init_tags:
        if tag in tagnames: emline_data[tag] = init_dict0[tag]

    # Mod on 12/02/2017 to include line_type
    # Mod on 02/03/2017 to include dl - spectral dispersion
    dict0 = {'data0':data0, 'emline_data':emline_data, 'fit_data0':fit_data0,
             'line_type':line_type, 'spec_file': os.path.basename(infile0),
             'x0': x0, 'dl': dl}

    # + on 10/02/2017
    if OH_file != None: dict0['OH_dict0'] = OH_dict0

    emline_data = fitting.main(dict0, out_pdf, out_fits, silent=silent,
                               verbose=verbose)

    if silent == False: log.info('### End read_data.main : '+systime())
#enddef

