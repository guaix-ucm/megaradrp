from __future__ import division
from __future__ import print_function

import argparse
import astropy.io.fits as fits
import numpy as np

import matplotlib.pyplot as plt
from numina.array.display.ximshow import ximshow
from numina.array.display.pause_debugplot import pause_debugplot
from numina.array.display.ximplotxy import ximplotxy
from numina.drps import get_system_drps


def cosinebell(n, fraction):
    """Return a cosine bell spanning n pixels, masking a fraction of pixels

    Parameters
    ----------
    n : int
        Number of pixels.
    fraction : float
        Length fraction over which the data will be masked.

    """

    mask = np.ones(n)
    nmasked = int(fraction * n)
    for i in range(nmasked):
        yval = 0.5 * (1 - np.cos(np.pi * float(i) / float(nmasked)))
        mask[i] = yval
        mask[n - i - 1] = yval

    return mask


def find_boxes(fitsfile, channels, nsearch, debugplot):
    """Refine boxes search around previous locations.

    Parameters
    ----------
    fitsfile : str
        FITS image where the new boxes will be measured.
    channels : tuple (integers)
        First and last channels (pixels in the x-direction) to extract
        median vertical cross section.
    nsearch : int
        Semi-width of the interval where each refined box location
        will be sought.
    debugplot : int
        Determines whether intermediate computations and/or plots
        are displayed:
        00 : no debug, no plots
        01 : no debug, plots without pauses
        02 : no debug, plots with pauses
        10 : debug, no plots
        11 : debug, plots without pauses
        12 : debug, plots with pauses
        21 : debug, additional plots without pauses
        22 : debug, additional plots with pauses

    Returns
    -------
    refined_boxes : numpy array
        Refined boxes locations.

    """

    # read the 2d image
    with fits.open(fitsfile) as hdulist:
        header = hdulist[0].header
        image2d = hdulist[0].data
    naxis2, naxis1 = image2d.shape
    print('>>> NAXIS1:', naxis1)
    print('>>> NAXIS2:', naxis2)

    if naxis1 != 4096 or naxis2 != 4112:
        raise ValueError("Unexpected (NAXIS1,NAXIS2) dimensions")

    # get previous boxes for current VPH, INSMODE and INSCONF
    vph = header['vph']
    insmode = header['insmode']
    insconf = header['insconf']
    print('>>> VPH:', vph)
    print('>>> INSMODE:', insmode)
    print('>>> INSCONF:', insconf)
    previous_boxes = get_previous_boxes(vph, insmode, insconf)

    if debugplot in (21, 22):
        ximshow(image2d, show=True,
                title='initial twilight image', debugplot=debugplot)

    # extract cross section
    nc1 = channels[0]
    nc2 = channels[1]
    ycut = np.median(image2d[:, nc1:(nc2+1)], axis=1)
    xcut = np.arange(naxis2) + 1
    if debugplot in (21, 22):
        ximplotxy(xcut, ycut, debugplot=debugplot,
                  xlabel='y axis', ylabel='number of counts',
                  title=fitsfile + " [" + str(nc1) + "," + str(nc2) + "]")

    # initial manipulation
    ycut -= np.median(ycut)  # subtract median
    ycut /= np.max(ycut)  # normalise values
    ycut *= -1  # invert signal to convert minima in maxima
    mask = cosinebell(n=ycut.size, fraction=0.10)
    ycut *= mask
    if debugplot in (21, 22):
        ximplotxy(xcut, ycut, debugplot=debugplot,
                  xlabel='y axis', ylabel='reversed scale')

    # Fourier filtering
    xf = np.fft.fftfreq(xcut.size)
    yf = np.fft.fftpack.fft(ycut)
    cut = (np.abs(xf) > 0.10)
    yf_trimmed = np.copy(yf)
    yf_trimmed[cut] = 0.0
    ycut_filt = np.fft.ifft(yf_trimmed).real
    if debugplot in (21, 22):
        ax = ximplotxy(xf, yf.real, plottype='semilog',
                       xlim=(0., 0.51), show=False,
                       label='original', linestyle='dotted')
        ax.plot(xf, yf_trimmed.real, label='trimmed')
        ax.legend()
        plt.show(block=False)
        plt.pause(0.001)
        pause_debugplot(debugplot)

    refined_boxes = np.zeros(previous_boxes.size, dtype=int)
    for ibox, box in enumerate(previous_boxes):
        iargmax = ycut_filt[box - nsearch:box + nsearch + 1].argmax()
        refined_boxes[ibox] = xcut[iargmax + box - nsearch]

    offsets = np.copy(refined_boxes)
    offsets -= previous_boxes
    print('>>> Offsets, new - old (pixels):', offsets)
    print('>>> New boxes:')
    nboxes = refined_boxes.size
    for i, dum in enumerate(refined_boxes):
        if i == nboxes - 1:
            print(dum)
        else:
            print(str(dum) + ',')

    if debugplot % 10 != 0:
        ax = ximplotxy(xcut, ycut_filt, show=False,
                       xlabel='y axis', ylabel='reversed scale')
        ax.vlines(previous_boxes, ymin=1.1, ymax=1.3, colors='magenta')
        ax.vlines(refined_boxes, ymin=1.4, ymax=1.6, colors='green')

        plt.show(block=False)
        plt.pause(0.001)
        pause_debugplot(debugplot)


def get_previous_boxes(vph, insmode, insconf):
    """Get previous boxes for the VPH, INSMODE and INSCONF.

    Using numina and megaradrp functionality.

    """

    d = get_system_drps()
    mydrp = d.drps['MEGARA']
    ic = mydrp.configurations[insconf]
    boxdict = ic.get('pseudoslit.boxes_positions',
                     **{'vph': vph, 'insmode': insmode})
    return np.array(boxdict['positions'])


def main(args=None):
    # parse command-line options
    parser = argparse.ArgumentParser(prog='find_boxes')
    # positional parameters
    parser.add_argument("fitsfile",
                        help="FITS image",
                        type=argparse.FileType('r'))
    parser.add_argument("--channels",
                        help="Channel region to extract cross section ",
                        default=(1990, 2010),
                        type=int, nargs=2)
    parser.add_argument("--nsearch",
                        help="Semi-width of the search window",
                        default=20, type=int)
    parser.add_argument("--debugplot",
                        help="integer indicating plotting/debugging" +
                             " (default=10)",
                        type=int, default=12,
                        choices=[0, 1, 2, 10, 11, 12, 21, 22])

    args = parser.parse_args(args=args)

    find_boxes(args.fitsfile.name, args.channels, args.nsearch,
               args.debugplot)


if __name__ == "__main__":

    main()
