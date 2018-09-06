import numpy as np
from scipy.interpolate import interp1d


from numpy.polynomial import Chebyshev as Ch


import pkg_resources

# load the sensitivity curves into an easy-to-use interpolation object.
sensfuncs_fname = pkg_resources.resource_filename("tresmerge", "data/sensfuncs.npy")
wls, sens = np.load(sensfuncs_fname)
# the sensfuncs is a (2,51,2998) shaped object, whereby the first is the wavelengths, and the second is the sensitivity.

# the number of echelle orders in TRES
norders = 51

# create order-by-order interpolation objects
interps = [interp1d(wls[i], sens[i], kind="linear", fill_value="extrapolate", assume_sorted=True) for i in range(norders)]



def optimize_calibration_static(wl_cal, fl_cal, sigma_cal, wl_fixed, fl_fixed, sigma_fixed, amp, l_f, order=1, mu_GP=1.0):
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
        amp (float): the GP amplitude
        l_f (float): the GP length
        order (int): the Chebyshev order to use
        mu_GP (optional): the mean of the GP to assume.

    Returns:
        (np.array, np.array): a tuple of two data products. The first is the ``fl_cal`` vector, now calibrated. The second is the array of the Chebyshev coefficients, in case one wants to re-evaluate the calibration polynomials.
    '''

    wl0 = np.min(wl)
    wl1 = np.max(wl)

    N_A = len(wl)
    A = np.empty((N_A, N_A), dtype=np.float64)

    N_B = len(wl_fixed)
    B = np.empty((N_B, N_B), dtype=np.float64)

    C = np.empty((N_A, N_B), dtype=np.float64)

    matrix_functions.fill_V11_f(A, wl_cal, amp, l_f)
    matrix_functions.fill_V11_f(B, wl_fixed, amp, l_f)
    matrix_functions.fill_V12_f(C, wl_cal, wl_fixed, amp, l_f)

    # Add in sigmas
    A[np.diag_indices_from(A)] += sigma_cal**2
    B[np.diag_indices_from(B)] += sigma_fixed**2


    # Get a clean set of the Chebyshev polynomials evaluated on the input wavelengths
    T = []
    for i in range(0, order + 1):
        coeff = [0 for j in range(i)] + [1]
        Chtemp = Ch(coeff, domain=[wl0, wl1])
        Ttemp = Chtemp(wl_cal)
        T += [Ttemp]

    T = np.array(T)

    D = fl_cal[:,np.newaxis] * T.T


    # Solve for the calibration coefficients c0, c1

    # Find B^{-1}, fl_prime, and C_prime
    try:
        B_cho = cho_factor(B)
    except np.linalg.linalg.LinAlgError:
        print("Failed to solve matrix inverse. Calibration not valid.")
        raise

    fl_prime = mu_GP + np.dot(C, cho_solve(B_cho, (fl_fixed.flatten() - mu_GP)))
    C_prime = A - np.dot(C, cho_solve(B_cho, C.T))

    # Find {C^\prime}^{-1}
    CP_cho = cho_factor(C_prime)

    # Invert the least squares problem
    left = np.dot(D.T, cho_solve(CP_cho, D))
    right = np.dot(D.T, cho_solve(CP_cho, fl_prime))

    left_cho = cho_factor(left)

    # the coefficents, X = [c0, c1]
    X = cho_solve(left_cho, right)

    # Apply the correction
    fl_cor = np.dot(D, X)

    return fl_cor, X


            # optimize the calibration of "tweak" with respect to all other orders
            fl_cor, X = covariance.optimize_calibration_ST1(lwl0, lwl1, lwl_tweak, fl_tweak, fl_remain[indsort], gp, A, C, order=config["order_cal"])

            # since fl_cor may have actually have fewer pixels than originally, we can't just
            # stuff the corrected fluxes directly back into the array.
            # Instead, we need to re-evaluate the line on all the wavelengths
            # (including those that may have been masked)
            # using the chebyshev coefficients, and apply this.

            # Here is where we need to make sure that wl0 and wl1 are the same.
            T = []
            for k in range(0, config["order_cal"] + 1):
                pass
                coeff = [0 for j in range(k)] + [1]
                Chtemp = Ch(coeff, domain=[lwl0, lwl1])
                Ttemp = Chtemp(lwl_tweak_unmasked)
                T += [Ttemp]

            T = np.array(T)

            Q = fl_tweak_unmasked[:,np.newaxis] * T.T

            # Apply the correction
            fl_cor = np.dot(Q, X)

            # replace this epoch with the corrected fluxes
            fl_out[i] = fl_cor[:,0]
