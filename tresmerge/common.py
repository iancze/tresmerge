import numpy as np
from scipy.interpolate import interp1d

from scipy.linalg import cho_factor, cho_solve
from numpy.polynomial import Chebyshev as Ch


# import pkg_resources
#
# # load the sensitivity curves into an easy-to-use interpolation object.
# sensfuncs_fname = pkg_resources.resource_filename("tresmerge", "data/sensfuncs.npy")
# wls, sens = np.load(sensfuncs_fname)
# the sensfuncs is a (2,51,2998) shaped object, whereby the first is the wavelengths, and the second is the sensitivity.

# the number of echelle orders in TRES
# norders = 51
#
# # create order-by-order interpolation objects
# interps = [interp1d(wls[i], sens[i], kind="linear", fill_value="extrapolate", assume_sorted=True) for i in range(norders)]



def solve_flux(wl, fl, sigma, fl_synth, order=3):
    '''
    Determine the calibration parameters for this epoch of observations. Assumes all wl, fl arrays are 1D, and that the relative velocities between all epochs are zero.

    Args:
        wl0 (float) : left wl point to evaluate the Chebyshev
        wl1 (float) : right wl point to evaluate the Chebyshev
        wl_cal (np.array) : the wavelengths of the epoch to calibrate
        fl_cal (np.array) : the fluxes of the epoch to calibrate
        sigma_cal (np.array): the sigmas of the epoch to calibrate
        wl_fixed (np.array) : the 1D (flattened) array of the reference wavelengths
        fl_fixed (np.array) : the 1D (flattened) array of the reference fluxes
        sigma_fixed (np.array) : the 1D (flattened) array of the reference sigmas
        order (int): the Chebyshev order to use
        mu_GP (optional): the mean of the GP to assume.

    Returns:
        (np.array, np.array): a tuple of two data products. The first is the ``fl_cal`` vector, now calibrated. The second is the array of the Chebyshev coefficients, in case one wants to re-evaluate the calibration polynomials.
    '''

    wl0 = np.min(wl)
    wl1 = np.max(wl)

    # Get a clean set of the Chebyshev polynomials evaluated on the input wavelengths
    T = []
    for i in range(0, order + 1):
        coeff = [0 for j in range(i)] + [1]
        Chtemp = Ch(coeff, domain=[wl0, wl1])
        Ttemp = Chtemp(wl)
        T += [Ttemp]

    T = np.array(T)

    # AKA, dcal
    D = fl[:,np.newaxis] * T.T

    # Inverse covariance matrix for data
    C_inv = np.diag(1/sigma**2)

    # Invert the least squares problem
    left = np.dot(D.T, np.dot(C_inv, D))
    right = np.dot(D.T, np.dot(C_inv, fl_synth))

    left_cho = cho_factor(left)

    # the coefficents, X = [c0, c1]
    X = cho_solve(left_cho, right)

    # Apply the correction to the data
    fl_cor = np.dot(D, X)

    return fl_cor, X
