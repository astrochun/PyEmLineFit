"""
fitting
=======

Python 2.7 codes to perform fitting to nebular emission lines with Gaussian
profiles
"""

import sys, os

from chun_codes import systime
from chun_codes import match_nosort
from chun_codes import chun_crossmatch

from os.path import exists
import commands
from astropy.io import ascii as asc
from astropy.io import fits

import numpy as np
import array, time, sets

import matplotlib.pyplot as plt
import glob

from astropy.table import Table
    
def main(dict0, out_pdf, silent=False, verbose=True):

    '''
    Provide explanation for function here.

    Parameters
    ----------
    dict0 : dictionary
      Dictionary containing:
        'data0' : Array of 1-D spectra
        'emline_data' : Emission-line fitting results astropy table
        'fit_data0' : astropy table containing emission lines to be fitted
        'x0' : numpy array of wavelength in Angstrom

    dict0 : dictionary
      Dictionary containing:

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True
	  
    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 9 September 2017
    '''
    
    if silent == False: print '### Begin fitting.main | '+systime()

    data0       = dict0['data0']
    emline_data = dict0['emline_data']
    fit_data0   = dict0['fit_data0']
    x0          = dict0['x0']
    zspec0      = emline_data['ZSPEC']

    n_spec = len(data0)

    for ll in xrange(n_spec):
        y0 = data0[ll]

        l_mark = np.where(y0 > 0)[0]
        if ll == 0: print l_mark
        if len(l_mark) > 0:
            emline_data['LMIN'][ll]  = np.min(x0[l_mark])
            emline_data['LMAX'][ll]  = np.max(x0[l_mark])
            emline_data['LMIN0'][ll] = np.min(x0[l_mark])/(1.0+zspec0[ll])
            emline_data['LMAX0'][ll] = np.max(x0[l_mark])/(1.0+zspec0[ll])

    if silent == False: print '### End fitting.main | '+systime()
    return emline_data
#enddef

