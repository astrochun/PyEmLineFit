"""
read_data
====

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
    
def main(infile0, zspec0=None, silent=False, verbose=True):

    '''
    Main function for reading in data

    Parameters
    ----------
    infile0 : string
      Filename for FITS data that contains 1-D spectra

    zspec0 : list or numpy array
      Array of spectroscopic redshifts for each 1-D spectrum

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
    
    if silent == False: print '### End read_data.main | '+systime()
#enddef

