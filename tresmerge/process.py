#!/usr/bin/env python

from tresmerge.common import interps

from astropy.io import ascii
from astropy.table import Table

from specutils.io import read_fits


def main():

    import argparse

    parser = argparse.ArgumentParser(description="Pseudo flux-calibrate and merge TRES echelle orders.")
    parser.add_argument("fits", help="Name of the FITS file containing the spectrum to merge.")
    parser.add_argument("--outfile", default="merged-spec.txt", help="Name of the output file to write the merged echelle spectrum to.")
    parser.add_argument("--clobber", action="store_true", help="Overwrite any existing output file with new data.")
    parser.add_argument("-t", "--trim", type=int, default=6, help="How many pixels to trim from the front of the file. Default is 6")
    args = parser.parse_args()


    # returns all of the echelle orders as a list
    raw_list = read_fits.read_fits_spectrum1d(args.fits)
    header = fits.getheader(args.fits)
    target = header["OBJECT"]
    date = header["DATE"]

    exptime = header.get("EXPTIME", 1)
    n_orders = 51
    n_pix = 2304 - args.trim

    # Create empty dataset to store pseudo-fluxcaled spectra
    wls = np.empty((n_orders, n_pix))
    fls = np.empty((n_orders, n_pix))
    sigmas = np.empty((n_orders, n_pix))

    # loop through the array, put things into an array of wl, fls
    # trimming and flux-caling along the way
    for order in range(51):

        # select the astropy spectrum corresponding to this order
        # and truncate according to trim
        raw = raw_list[order]
        wl = raw.wavelength.value[args.trim:]
        fl = raw.flux.value[args.trim:]

        # calculate the uncertainty on this pixel as simply the sqrt of the number of counts
        sigma = np.sqrt(fl)

        # select the sensitivity curve interpolation object specific to this order
        interp = interps[order]

        # calculate the instrument sensitivity at these wavelengths
        Sfunc = interp(wl)

        flcal = fl / EXPTIME * 10**(-0.4 * Sfunc)
        sigmacal = sigma / EXPTIME * 10**(-0.4 * (Sfunc)

        wls[order] = wl
        fls[order] = flcal
        sigmas[order] = sigmacal

    # make some plots to examine what the spectrum looks like at the overlap regions
    with PdfPages('output.pdf', metadata={'creator':'tresmerge', 'author':'Ian Czekala', 'title':"{:} - {:}".format(target, date)}) as pdf:

        for order in range(n_orders):

            order_min = order - 1 if order > 0 else 0
            order_max = order + 1 if order < (n_orders - 2) else (norders - 1)

            wl_left = wls[order_min]
            fl_left = fls[order_min]

            wl_center = wls[order]
            fl_center = fls[order]

            wl_right = wls[order_max]
            fl_right = fls[order_max]

            # As many times as you like, create a figure fig and save it:
            fig, ax = plt.subplots(nrows=1, figsize=(8,4))

            ax.plot(wl_left, fl_left, color="0.6")
            ax.plot(wl_right, fl_right, color="0.6")
            ax.plot(wl_center, fl_center, color="k")

            # save figure to large PDF file
            pdf.savefig(fig)

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
