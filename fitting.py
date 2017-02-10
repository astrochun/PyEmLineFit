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

from matplotlib.backends.backend_pdf import PdfPages # + on 10/02/2017

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
    Created by Chun Ly, 9 February 2017
    Modified by Chun Ly, 10 February 2017
     - Write PDF to out_pdf and show panels for each emission line
    '''
    
    if silent == False: print '### Begin fitting.main | '+systime()

    data0       = dict0['data0']
    emline_data = dict0['emline_data']
    fit_data0   = dict0['fit_data0']
    x0          = dict0['x0']
    zspec0      = emline_data['ZSPEC']

    n_spec = len(data0)

    pp = PdfPages(out_pdf) # + on 10/02/2017

    for ll in xrange(n_spec):
        if verbose == True:
            print '### Working on line=%04i zspec : %5.3f' % (ll, zspec0[ll])

        y0 = data0[ll]

        # Fill in LMIN/LMAX values, range of spectrum
        l_mark = np.where(y0 > 0)[0]
        if len(l_mark) > 0:
            lmin, lmax = np.min(x0[l_mark]), np.max(x0[l_mark])
            emline_data['LMIN'][ll]  = lmin
            emline_data['LMAX'][ll]  = lmax
            emline_data['LMIN0'][ll] = np.min(x0[l_mark])/(1.0+zspec0[ll])
            emline_data['LMAX0'][ll] = np.max(x0[l_mark])/(1.0+zspec0[ll])

        # + on 10/02/2017
        z_lines = (1.0+zspec0[ll]) * fit_data0['lines0'] # redshifted lines
        x_min   = (1.0+zspec0[ll]) * fit_data0['x_min0']
        xwidth  = (1.0+zspec0[ll]) * fit_data0['xwidth0']

        # + on 10/02/2017
        mask_lines_min = z_lines - 10.
        mask_lines_max = z_lines + 10.

        # Array of wavelengths affected by emission lines | + on 10/02/107
        line_flag = np.zeros(len(y0))
        for kk in xrange(len(fit_data0)):
            mark2 = np.where((x0 >= mask_lines_min[kk]) &
                             (x0 <= mask_lines_max[kk]))[0]
            if len(mark2) > 0: line_flag[mark2] = 1

        in_spec = np.where((z_lines >= lmin) & (z_lines <= lmax))[0]

        # + on 10/02/2017
        if len(in_spec) == 0:
            print '### No lines available'
        else:
            fig0, ax_arr = plt.subplots(nrows=4, ncols=2)
            for ss in xrange(len(in_spec)):
                panel_check = 1 if ss == 0 else 0
                if ss != 0:
                    panel_check = fit_data0['panels'][in_spec[ss]] - \
                                  fit_data0['panels'][in_spec[ss-1]]

                t_col = fit_data0['panels'][in_spec[ss]] % 2
                t_row = fit_data0['panels'][in_spec[ss]] / 2
                if panel_check:
                    t_ax0 = ax_arr[t_row,t_col]
                    t_ax0.plot(x0,y0)
                    xra = [x_min[in_spec[ss]],
                           x_min[in_spec[ss]]+2*xwidth[in_spec[ss]]]
                    t_ax0.set_xlim(xra)
            #endfor
            fig0.savefig(pp, format='pdf')
    #endfor

    pp.close()

    if silent == False: print '### End fitting.main | '+systime()
    return emline_data
#enddef

