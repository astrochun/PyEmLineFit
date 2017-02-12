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
import matplotlib.patches as patches # + on 10/02/2017

import glob

from astropy.table import Table

from pylab import subplots_adjust # + on 10/02/2017

from matplotlib.backends.backend_pdf import PdfPages # + on 10/02/2017

def draw_OH(OH_dict0, t_ax0, xra, yra, silent=True, verbose=False):
    '''
    Shade wavelengths affected by OH night skylines

    Parameters
    ----------
    OH_dict0 : dict
      Dictionary containing:
        'data0' : Array of 1-D spectra
        'emline_data' : Emission-line fitting results astropy table
        'fit_data0' : astropy table containing emission lines to be fitted
        'x0' : numpy array of wavelength in Angstrom

    t_ax0 : matplotlib.axes
      Plotting axes from plt.subplots

    xra : array or list
      Limit of x-axis for plotting

    yra : array or list
      Limit of y-axis for plotting

    silent : boolean
      Turns off stdout messages. Default: True

    verbose : boolean
      Turns on additional stdout messages. Default: False

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 10 February 2017
    Modified by Chun Ly, 11 February 2017
     - Minor bug in patches call with alpha and color
    '''

    OH_xmin0 = OH_dict0['OH_xmin0']
    OH_xmax0 = OH_dict0['OH_xmax0']

    OH_in_reg = np.where((OH_xmin0 > xra[0]) & (OH_xmin0 < xra[1]))[0]
    if len(OH_in_reg) > 0:
        for oo in xrange(len(OH_in_reg)):
            o_idx = OH_in_reg[oo]
            t_max = np.min([OH_xmax0[o_idx], xra[1]])
            t_min = OH_xmin0[o_idx]
            dx    = t_max - t_min
            dy    = yra[1] - yra[0]
            # Minor bug fixed on 11/02/2017
            t_ax0.add_patch(patches.Rectangle((t_min, yra[0]), dx, dy,
                                              alpha=0.5, facecolor='k',
                                              edgecolor='none'))
#enddef

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

    out_pdf : string
      Filename for output PDF file

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
     - Adjust plotting to include label on top, draw emission lines,
       draw zero
     - Add call to draw_OH()
     - Added use of [line] for indexing of 1-D spectra to fit
    Modified by Chun Ly, 12 February 2017
     - Handle different input types ('SLIT', 'AP') for title labeling
     - Annotate emission lines for each panel in the upper right
     - Define [no_spec] to properly handle not plotting certain subplots panel
    '''
    
    if silent == False: print '### Begin fitting.main | '+systime()

    data0       = dict0['data0']
    emline_data = dict0['emline_data']
    fit_data0   = dict0['fit_data0']
    spec_file   = dict0['spec_file'] # + on 10/02/2017
    x0          = dict0['x0']
    zspec0      = emline_data['ZSPEC']
    line        = emline_data['LINE'] # + on 11/02/2017

    # Get OH skyline info | + on 10/02/2017
    has_OH = 0
    if 'OH_dict0' in dict0.keys():
        OH_flag0 = dict0['OH_dict0']['OH_flag0']
        OH_xmin0 = dict0['OH_dict0']['OH_xmin0']
        OH_xmax0 = dict0['OH_dict0']['OH_xmax0']
        has_OH = 1
    else:
        OH_flag0 = np.zeros(x0, dtype=np.int8)

    n_spec = len(data0)
    n_line = len(line) # + on 11/02/2017

    if silent == False: print '### Output PDF file : ', out_pdf
    pp = PdfPages(out_pdf) # + on 10/02/2017
    nrows, ncols = 4, 2

    for ll in xrange(10): #n_line):
        if verbose == True:
            print '### Working on line=%04i zspec : %.4f' % (line[ll],
                                                             zspec0[ll])

        y0 = data0[line[ll]-1]

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

        # Emission lines not in spectral coverage | + on 12/02/2017
        no_spec = np.where((z_lines < lmin) | (z_lines > lmax))[0]

        # + on 10/02/2017
        if len(in_spec) == 0:
            print '### No lines available'
        else:
            fig0, ax_arr = plt.subplots(nrows=nrows, ncols=ncols)

            # + on 10/02/2017, Mod on 12/02/2017
            label0 = spec_file
            if 'SLIT' in emline_data.colnames:
                label0 += ' SLIT=%s' % emline_data['SLIT'][ll]
            if 'AP' in emline_data.colnames:
                label0 += ' AP=%s' % emline_data['AP'][ll]
            label0 += ' line=%04i' % line[ll]
            label0 += ' z=%6.4f %.1f-%.1f' % (zspec0[ll], lmin, lmax)
            ax_arr[0,0].set_title(label0, loc=u'left', fontsize=14)

            for ss in xrange(len(in_spec)):
                s_idx = in_spec[ss] # + on 10/02/2017
                panel_check = 1 if ss == 0 else 0

                ss_panel = fit_data0['panels'][s_idx] # + on 12/02/2017
                if ss != 0:
                    panel_check = ss_panel - fit_data0['panels'][in_spec[ss-1]]

                t_col, t_row = ss_panel % ncols, ss_panel / ncols

                t_ax0 = ax_arr[t_row,t_col] #Moved up on 10/02/2017
                if panel_check:
                    xra = [x_min[s_idx], x_min[s_idx]+2*xwidth[s_idx]]

                in_range  = np.where((x0 >= xra[0]) & (x0 <= xra[1]))[0]
                in_range2 = np.where((x0 >= xra[0]) & (x0 <= xra[1]) &
                                     (OH_flag0 == 0))[0]

                # Mod on 10/02/2017
                if panel_check: #Do this once for each panel
                    t_ax0.plot(x0, y0, 'k', linewidth=0.75, zorder=5)
                    t_ax0.set_xlim(xra)

                    # Adjust y-range
                    t1  = np.min(y0[in_range2])
                    yra = [1.1*t1 if t1 < 0 else 0.0, 1.1*np.max(y0[in_range2])]
                    t_ax0.set_ylim(yra)
                    t_ax0.yaxis.set_ticklabels([])
                    t_ax0.tick_params(axis='both', which='major', labelsize=10)
                    t_ax0.minorticks_on()

                    # Draw OH skyline | + on 10/02/2017
                    if has_OH: draw_OH(dict0['OH_dict0'], t_ax0, xra, yra)

                    # Draw horizontal line at y=0
                    t_ax0.axhline(y=0.0, linewidth=1, color='b', zorder=1)

                    # Annotate emission lines in the upper right panel | + on 12/02/2017
                    in_panel = np.where(fit_data0['panels'][in_spec] == ss_panel)[0]
                    in_panel = in_spec[in_panel]
                    line_annot = '\n'.join([a for a in fit_data0['lines_greek'][in_panel]])
                    bbox_props = dict(boxstyle="square,pad=0.3", fc="white", alpha=0.9,
                                      ec="none")
                    t_ax0.annotate(line_annot, (0.95,0.95), xycoords='axes fraction',
                                   ha='right', va='top', bbox=bbox_props, zorder=6)

                # Draw vertical lines for axes | + on 10/02/2017
                t_ax0.axvline(x=z_lines[s_idx], linewidth=1, color='b',
                              zorder=1)
            #endfor

            # Remove plotting of non-use panels | + on 10/02/2017
            if len(no_spec) > 0:
                no_panel = list(set(fit_data0['panels'][no_spec])) # Remove duplicates
                for npl in no_panel:
                    t_col, t_row = npl % ncols, npl / ncols
                    ax_arr[t_row,t_col].axis('off')

            # + on 10/02/2017
            subplots_adjust(left=0.025, bottom=0.025, top=0.975, right=0.975,
                            wspace=0.05, hspace=0.10)
            fig0.set_size_inches(8, 11)
            fig0.savefig(pp, format='pdf')
        #endelse
    #endfor

    pp.close()

    if silent == False: print '### End fitting.main | '+systime()
    return emline_data
#enddef

