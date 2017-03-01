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

import read_data
import fitting
import execute
import test

reload(read_data)
reload(fitting)
reload(execute)
reload(test)

#__all__ = ['read_data']
