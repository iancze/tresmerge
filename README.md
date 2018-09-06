# tresmerge
Pseudo-flux calibrate and merge TRES echelle spectra

Requires Python 3.

Required packages:

* numpy
* scipy
* matplotlib
* astropy
* astropy/specutils: https://github.com/astropy/specutils

## Installation

Download this package, `cd` to the directory, and run

  $python setup.py install

or

  $python3 setup.py install

Make sure that `python` or `python3` points to the Python 3 version you want to use. You can check that this is the correct version by first entering the Python interpreter, and you should see something like

  $ python
  Python 3.6.3 |Anaconda custom (64-bit)| (default, Nov  3 2017, 19:19:16)
  [GCC 7.2.0] on linux
  Type "help", "copyright", "credits" or "license" for more information.
  >>>

These scripts have not been tested with Python 2.x and will likely not work. It's recommended you upgrade to Python 3.
