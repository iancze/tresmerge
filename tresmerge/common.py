import numpy as np
from scipy.interpolation import interp1d
import pkg_resources

# load the sensitivity curves into an easy-to-use interpolation object.
sensfuncs_fname = pkg_resources.resource_filename("tresmerge", "data/sensfuncs.npy")
wls, sens = np.load(sensfuncs_fname)
# the sensfuncs is a (2,51,2998) shaped object, whereby the first is the wavelengths, and the second is the sensitivity.

# the number of echelle orders in TRES
norders = 51

# create order-by-order interpolation objects
interps = [interp1d(wls[i], sens[i], kind="linear", fill_value="extrapolate", assume_sorted=True) for i in range(norders)]
