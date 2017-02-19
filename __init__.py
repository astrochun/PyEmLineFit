"""
PyEmLineFit
===========

Package of codes created to fit nebular emission lines in 1-D spectra

Requirements:
 pyspeckit (https://pyspeckit.bitbucket.io/?)

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
