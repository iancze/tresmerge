#!/usr/bin/env python


import numpy as np
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from astropy.io import ascii, fits
from astropy.table import Table

from specutils.io import read_fits

c_ang = 2.99792458e18 #A s^-1
c_kms = 2.99792458e5 #km s^-1

n_air = 1.000277
c_ang_air = c_ang/n_air
c_kms_air = c_kms/n_air

def main():

    import argparse

    parser = argparse.ArgumentParser(description="Merge TRES echelle orders.")
    parser.add_argument("rfits", help="Name of the FITS file containing the RAW spectrum to merge.")
    parser.add_argument("bfits", help="Name of the FITS file containing the BLAZE-corrected spectrum to merge.")
    parser.add_argument("template", help="Name of the FITS file containing the Kurucz template spectrum.")
    parser.add_argument("--outfile", default="merged-spec.txt", help="Name of the output file to write the merged echelle spectrum to.")
    parser.add_argument("--clobber", action="store_true", help="Overwrite any existing output file with new data.")
    parser.add_argument("-t", "--trim", type=int, default=6, help="How many pixels to trim from the front of the file. Default is 6")
    parser.add_argument("--shift", default=0.0, type=float, help="Doppler shift the synthetic spectrum by this amount (in km/s) before doing the merge process. This may help if the target star has an exceptionally high radial velocity. Positive velocities correspond to redshifting the template.")
    parser.add_argument("--poly-order", default=4, type=int, help="The order Chebyshev polynomial used to flatten each echelle order. 0 = constant, 1 = line, 2 = parabola, 3 = cubic, ... etc.")
    args = parser.parse_args()

    ######################
    # TEMPLATE PROCESSING
    # read the template
    hdul_t = fits.open(args.template)
    fl_t = hdul_t[0].data
    hdr_t = hdul_t[0].header

    # calculate the wavelength vector from the quantities in the FITS header
    num = len(fl_t)
    p = np.arange(num)
    w1 = hdr_t['CRVAL1']
    dw = hdr_t['CDELT1']
    wl_t = 10 ** (w1 + dw * p)
    hdul_t.close()

    dopp_shift = np.sqrt((c_kms_air + args.shift) / (c_kms_air - args.shift))
    wl_t = wl_t * dopp_shift

    #Also, we should convert from f_nu to f_lam
    fl_t *= c_ang / wl_t**2 #Convert from f_nu to f_lambda
    fl_t /= np.average(fl_t) #divide by the mean flux, so avg(f) = 1

    # create an interpolation object to the synthetic spectrum
    interp = interp1d(wl_t, fl_t, kind="cubic", bounds_error=True)
    ######################


    ######################
    # READING IN THE DATA
    ######################
    # returns all of the echelle orders as a list
    raw_list = read_fits.read_fits_spectrum1d(args.rfits)
    blaze_list = read_fits.read_fits_spectrum1d(args.bfits)
    header = fits.getheader(args.rfits)
    target = header["OBJECT"]
    date = header["DATE"]

    n_orders = 51
    n_pix = 2304 - args.trim

    # Create empty dataset to store pseudo-fluxcaled spectra
    wls = np.empty((n_orders, n_pix))
    fls = np.empty((n_orders, n_pix))
    fls_t = np.empty((n_orders, n_pix)) # to store the interpolated model fluxes
    sigmas = np.empty((n_orders, n_pix))

    # loop through the array, put things into an array of wl, fls
    # trimming and flux-caling along the way
    for order in range(51):
        # select the astropy spectrum corresponding to this order
        # and truncate according to trim
        raw = raw_list[order]
        blaze = blaze_list[order]

        wl = raw.wavelength.value[args.trim:]
        fl_r = raw.flux.value[args.trim:]
        # Where the ratio values are 0, just set it to 1, since the noise will be 0 here too.
        fl_r[fl_r==0] = 1.0
        fl_b = blaze.flux.value[args.trim:]

        ratio = fl_b / fl_r

        # calculate the uncertainty on this pixel as simply the sqrt of the number of counts,
        # but in units of the blaze function
        sigma = ratio * np.sqrt(fl_r)

        wls[order] = wl
        fls[order] = fl_b
        sigmas[order] = sigma

        try:
            fls_t[order] = interp(wl)
        except ValueError:
            print("Cannot interpolate synthetic spectrum for echelle order {:}. Assuming flat continuum for now, but you may want to discard the spectrum in this order.".format(d))
            fls_t[order] = np.ones_like(wl)


    ############################################
    # FLATTENING THE DATA TO MATCH THE TEMPLATE
    ###########################################

    # Norm the average level of all orders to match the level of the synthetic spectrum
    fl_median = np.median(fls, axis=1)
    flt_median = np.median(fls_t, axis=1)
    factors = flt_median / fl_median
    # update these to be normalized to the synthetic spectra
    fls = fls * factors[:,np.newaxis]
    sigmas = sigmas * factors[:,np.newaxis]

    # get the coefficients for the polynomial which makes it match the synthetic spectrum.

    ######################
    # PLOTTING THE DATA
    ######################
    # make some plots to examine what the spectrum looks like at the overlap regions
    with PdfPages('output.pdf', metadata={'Creator':'tresmerge', 'Author':'Ian Czekala', 'Title':"{:} - {:}".format(target, date)}) as pdf:

        for order in range(n_orders):

            order_min = order - 1 if order > 0 else 0
            order_max = order + 1 if order < (n_orders - 2) else (n_orders - 1)

            wl_left = wls[order_min]
            fl_left = fls[order_min]

            wl_center = wls[order]
            fl_center = fls[order]
            fl_t_center = fls_t[order]

            wl_right = wls[order_max]
            fl_right = fls[order_max]

            # As many times as you like, create a figure fig and save it:
            fig, ax = plt.subplots(nrows=1, figsize=(8,4))

            ax.plot(wl_left, fl_left, color="0.6")
            ax.plot(wl_right, fl_right, color="0.6")
            ax.plot(wl_center, fl_center, color="k")
            ax.set_xlabel(r"$\lambda$\,[\AA]")
            ax.set_ylabel(r"$f_\lambda\,[\propto$\,counts]")
            ax.set_title("Order {:d}".format(order + 1))

            ax.plot(wl_center, fl_t_center)

            # save figure to large PDF file
            pdf.savefig(fig)
            plt.close('all')
    #
    #
    # # construct an astropy table with the merged files
    # t = Table(, names=["wl", "fl", "sigma"])
    #
    # # put the header into the ECSV header
    # t.meta = header
    #
    # # specify how much accuracy to use for each column
    # formats = {"wl":"%.3f", "fl":"%.3e", "sigma":"%.3e"}
    #
    # # write the result to file
    # ascii.write(t, arg.outfile, overwrite=args.clobber, format="ecsv", formats=formats)

if __name__=="__main__":
    main()
