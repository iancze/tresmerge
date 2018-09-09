#!/usr/bin/env python

import datetime
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline, UnivariateSpline

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from astropy.io import ascii, fits
from astropy.table import Table

from specutils.io import read_fits

from tresmerge.common import solve_flux

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
    parser.add_argument("--plot", action="store_true", help="Make a set of plots of the merged spectra.")
    parser.add_argument("-t", "--trim", type=int, default=6, help="How many pixels to trim from the front of the file. Default is 6")
    parser.add_argument("--shift", default=0.0, type=float, help="Doppler shift the synthetic spectrum by this amount (in km/s) before doing the merge process. This may help if the target star has an exceptionally high radial velocity. Positive velocities correspond to redshifting the template.")
    parser.add_argument("--poly-order", default=3, type=int, help="The order Chebyshev polynomial used to flatten each echelle order. 0 = constant, 1 = line, 2 = parabola, 3 = cubic, ... etc.")
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
    interp = InterpolatedUnivariateSpline(wl_t, fl_t, k=2, ext=2)
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
    fls_cor = np.empty((n_orders, n_pix))
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
        fl_b[fl_b==0] = 1.0

        ratio = fl_b / fl_r

        # calculate the uncertainty on this pixel as simply the sqrt of the number of counts,
        # but in units of the blaze function
        sigma = ratio * np.sqrt(fl_r + 0.01 * np.median(fl_r)) # to avoid very low vals

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

    for order in range(n_orders):
        # get the coefficients for the polynomial which makes it match the synthetic spectrum.
        fl_cor, X = solve_flux(wls[order], fls[order], sigmas[order], fls_t[order], order=args.poly_order)
        fls_cor[order,:] = fl_cor

    ######################
    # Smoothing the overlap
    ######################

    # start with an existing stub, and add onto it with each iteration
    wl_all = wls[0]
    fl_all = fls_cor[0]
    sigma_all = sigmas[0]

    for order in range(1, n_orders):

        # Find all wavelengths in the stub that have a value larger than the shortest wavelength in the red order
        l_ind = (wl_all > np.min(wls[order]))

        # Find all wavelengths in the red order that have a value smaller than the largest wavelength in the stub
        r_ind = (wls[order] < np.max(wl_all))

        # if either of these is less than 2, then there isn't enough overlap to actually merge and we should
        # just take add on the new order to the stub without doing any fancy smoothing
        if (np.sum(l_ind) < 2) or (np.sum(r_ind) < 2):

            wl_all = np.concatenate((wl_all, wls[order]))
            fl_all = np.concatenate((fl_all, fls_cor[order]))
            sigma_all = np.concatenate((sigma_all, sigmas[order]))

        else:

            wl_smooth = wl_all[l_ind]

            # Grab both wavelengths, fluxes, and sigmas in the overlap region
            wl_comb = np.concatenate([wl_smooth, wls[order][r_ind]])
            fl_comb = np.concatenate([fl_all[l_ind], fls_cor[order][r_ind]])
            sigma_comb = np.concatenate([sigma_all[l_ind], sigmas[order][r_ind]])

            # sort according to wl
            indsort = np.argsort(wl_comb)
            wl_sorted = wl_comb[indsort]
            fl_sorted = fl_comb[indsort]
            sigma_sorted = sigma_comb[indsort]

            # construct a spline
            spline = UnivariateSpline(wl_sorted, fl_sorted, w=1/sigma_sorted)

            # choose the left_wl as the wl to interpolate onto
            fl_smooth = spline(wl_smooth)

            # transfer the sigmas, too, as the sqrt of sum of squares of the adjacent pixels
            # go through each wl in wl_all[lind], find it's arguement in wl_sorted, and average the two
            # nearest pixels in sigma_sorted
            n_comb = len(sigma_sorted)
            n_smooth = len(wl_smooth)

            sigma_smooth = np.empty(n_smooth)

            for i in range(n_smooth):
                w = wl_smooth[i]

                # find where this appears in the combined array
                arg = np.searchsorted(wl_sorted, w)

                # choose this is the first value as long as it isn't the rightmost value
                if arg < (n_smooth - 1):
                    sigma_left = arg
                    sigma_right = arg + 1
                else:
                    sigma_left = arg - 1
                    sigma_right = arg

                sigma_smooth[i] = np.sqrt(sigma_left**2 + sigma_right**2)

            # append all of the pieces together to grow the merged spectrum until we're done adding orders
            wl_all = np.concatenate((wl_all, wls[order][~r_ind]))
            fl_all = np.concatenate((fl_all[~l_ind], fl_smooth, fls_cor[order][~r_ind]))
            sigma_all = np.concatenate((sigma_all[~l_ind], sigma_smooth, sigmas[order][~r_ind]))


    if args.plot:
        ######################
        # PLOTTING THE DATA
        ######################
        # make some plots to examine what the spectrum looks like at the overlap regions
        with PdfPages('output.pdf') as pdf:

            for order in range(n_orders):

                order_min = order - 1 if order > 0 else 0
                order_max = order + 1 if order < (n_orders - 2) else (n_orders - 1)

                wl_left = wls[order_min]
                fl_left = fls[order_min]
                fl_t_left = fls_t[order_min]
                fl_cor_left = fls_cor[order_min]

                wl_center = wls[order]
                fl_center = fls[order]
                fl_t_center = fls_t[order]
                fl_cor_center = fls_cor[order]

                wl_right = wls[order_max]
                fl_right = fls[order_max]
                fl_t_right = fls_t[order_max]
                fl_cor_right = fls_cor[order_max]

                # As many times as you like, create a figure fig and save it:
                fig, ax = plt.subplots(nrows=3, figsize=(8.5,12), sharex=True)

                ax[0].plot(wl_left, fl_left, color="0.6")
                ax[0].plot(wl_right, fl_right, color="0.6")
                ax[0].plot(wl_center, fl_center, color="k", label="data")
                ax[0].set_title("Order {:d}".format(order + 1))

                ax[0].plot(wl_center, fl_cor_center, "r", label="cor. data")
                ax[0].plot(wl_center, fl_t_center, "b", label="template")
                ax[0].legend(loc="best", fontsize="x-small")

                # The polynomial-corrected spectra, with the synthetic model on top

                ax[1].plot(wl_left, fl_cor_left, "r", label="cor. data edge")
                ax[1].plot(wl_right, fl_cor_right, "r")
                ax[1].plot(wl_center, fl_cor_center, "g", alpha=0.7, label="cor. data center")
                ax[1].legend(loc="best", fontsize="x-small")

                ax[1].plot(wl_left, fl_t_left, "b", label="template")
                ax[1].plot(wl_center, fl_t_center, "b")
                ax[1].plot(wl_right, fl_t_right, "b")

                # ax[1].plot(wl_t, fl_t, "k", lw=0.5)
                ax[1].set_title("Flattened orders")


                # The final merged spectrum
                ax[2].plot(wl_all, fl_all, "b", label="merged data")
                ax[2].legend(loc="best", fontsize="x-small")
                ax[2].set_xlim(np.min(wl_left), np.max(wl_right))
                ax[2].set_ylim(ax[1].get_ylim())
                ax[2].set_title("Merged spectrum")


                ax[0].set_ylabel(r"$f_\lambda\,[\propto$\,counts]")
                ax[1].set_ylabel(r"$f_\lambda\,[\propto$\,counts]")
                ax[2].set_ylabel(r"$f_\lambda\,[\propto$\,counts]")

                ax[2].set_xlabel(r"$\lambda$\,[\AA]")
                # save figure to large PDF file
                pdf.savefig(fig)

                plt.close('all')


            # We can also set the file's metadata via the PdfPages object:
            d = pdf.infodict()
            d['Title'] = "{:} - {:}".format(target, date)
            d['Author'] = "Ian Czekala"
            d['CreationDate'] = datetime.datetime(2009, 11, 13)
            d['ModDate'] = datetime.datetime.today()

    # sort all of the wl, fl, and sigmas just for good measure
    indsort = np.argsort(wl_all)
    wl_all = wl_all[indsort]
    fl_all = fl_all[indsort]
    sigma_all = sigma_all[indsort]

    # construct an astropy table with the merged files
    t = Table((wl_all, fl_all, sigma_all), names=["wl", "fl", "sigma"])

    for (key,value) in header.items():
        t.meta[key] = value

    # specify how much accuracy to use for each column
    formats = {"wl":"%.3f", "fl":"%.3e", "sigma":"%.3e"}

    # write the result to file
    ascii.write(t, args.outfile, overwrite=args.clobber, format="ecsv", formats=formats)

if __name__=="__main__":
    main()
