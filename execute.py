"""
execute
====

Provide description for code here.
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

import collections

import read_data

def deep2(silent=False, verbose=True):

    '''
    Execute PyEmLinefit for the DEEP2 DR4 data

    Parameters
    ----------

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True
	  
    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 11 February 2017
    Modified by Chun Ly, 12 February 2017
     - Change order of init_dict0 so dtype is first, following format
    '''
    
    if silent == False: print '### Begin execute.deep2 | '+systime()

    deep2_path = '/Users/cly/data/DEEP2/DR4/'

    infiles   = glob.glob(deep2_path+'DEEP2_2D_Field?.f3.fits')
    cat_files = glob.glob(deep2_path+'DEEP2_Field?_all_line_fit.cat.fits')

    OH_files  = glob.glob(deep2_path+'DEEP2_Field?.OH.txt')

    pdf_files = [str0.replace('.fits','.ELF.pdf') for str0 in infiles]

    for ff in xrange(len(infiles)):
        cat_data = fits.getdata(cat_files[ff])

        objno = cat_data.OBJNO
        zspec = cat_data.ZBEST
        line  = 1+np.arange(len(cat_data))
        slit  = ['%4i.%03i' % (a,b) for a,b in
                 zip(cat_data.MASK, cat_data.SLIT)]

        dtype0 = ['i', 'S8', 'i', 'f8']
        init_dict0 = collections.OrderedDict()
        init_dict0['dtype'] = dtype0
        init_dict0['OBJNO'] = objno
        init_dict0['SLIT']  = slit
        init_dict0['LINE']  = line
        init_dict0['ZSPEC'] = zspec
        
        read_data.main(infiles[ff], init_dict0, OH_file=OH_files[ff],
                       resol='high', out_pdf=pdf_files[ff])
    if silent == False: print '### End execute.deep2 | '+systime()
#enddef

