"""
PyEmLineFit
===========

Package of codes created to fit nebular emission lines in 1-D spectra

Requirements:
 astropy    >=v1.3
 chun_codes
 matplotlib >= v1.5.3
 numpy      >= v1.11.1
 pyspeckit  >= v0.1.19 (https://pyspeckit.bitbucket.io/?)

How to use the documentation
----------------------------
import PyEmLineFit as ELF
"""

# + on 20-21/03/2017
cols = ['LAMBDA', 'ZSPEC', 'PEAK', 'SIGMA', 'Y0', 'FLUX', 'FLUX_ERR',
        'FLUX_DATA', 'NOISE', 'SNR']
c_dtype = ['f8','f8','e','f8','e','e','e','e','e','f8']

__all__ = ['cols', 'c_dtype']

__version__ = '0.1'

import read_data
import fitting
import execute
import test

reload(read_data)
reload(fitting)
reload(execute)
reload(test)

