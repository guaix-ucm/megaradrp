#
# Copyright 2011-2015 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# Megara DRP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Megara DRP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Megara DRP.  If not, see <http://www.gnu.org/licenses/>.
#

from __future__ import print_function

from astropy.io import fits
import numpy as np

from numina.array.trace.extract import extract_simple_rss
from numina.array.utils import wc_to_pix_1d as wcs_to_pix_1d

# row / column
_binning = {'11': [1, 1], '21': [1, 2], '12': [2, 1], '22': [2, 2]}
_direc = ['normal', 'mirror']


def trim_and_o(image, out='trimmed.fits', direction='normal', bins='11'):
    """Trim a MEGARA image with overscan."""

    with fits.open(image) as hdul:
        hdu = trim_and_o_hdu(hdul[0], direction, bins)
        hdu.writeto(out, clobber=True)


def trim_and_o_hdu(hdu, direction='normal', bins = '11'):
    """Trim a MEGARA HDU with overscan."""

    finaldata = trim_and_o_array(hdu.data, direction=direction, bins=bins)

    hdu.data = finaldata
    return hdu


def trim_and_o_array(array, direction='normal', bins='11'):
    """Trim a MEGARA array with overscan."""

    if direction not in _direc:
        raise ValueError("%s must be either 'normal' or 'mirror'" % direction)

    if direction == 'normal':
        direcfun = lambda x: x
    else:
        direcfun = np.fliplr

    if bins not in _binning:
        raise ValueError("%s must be one if '11', '12', '21, '22'" % bins)

    OSCANW = 100
    PSCANW = 50
    H_X_DIM = 2048
    H_Y_DIM = 2056

    bng = _binning[bins]

    nr2 = H_Y_DIM * 2 / bng[0]
    nc2 = H_X_DIM * 2 / bng[1]

    nr = H_Y_DIM / bng[0]

    oscan2 = OSCANW / bng[0]
    psc1 = PSCANW / bng[0]

    finaldata = np.empty((nr2, nc2), dtype='float32')
    finaldata[:nr, :] = direcfun(array[:nr, psc1:nc2 + psc1])
    finaldata[nr:, :] = direcfun(array[nr + oscan2:, psc1:nc2 + psc1])
    return finaldata


def apextract(data, trace):
    """Extract apertures."""
    rss = np.empty((trace.shape[0], data.shape[1]), dtype='float32')
    for idx, r in enumerate(trace):
        l = r[0]
        r = r[2] + 1
        sl = (slice(l, r), )
        m = data[sl].sum(axis=0)
        rss[idx] = m
    return rss


def apextract_tracemap(data, tracemap):
    """Extract apertures using a tracemap."""
    
    # FIXME: a little hackish
    
    pols = [np.poly1d(t['fitparms']) for t in tracemap]

    borders = []

    # These are polynomial, they can be summed

    # Estimate left border for the first trace
    pix_12 = 0.5 * (pols[0] + pols[1])
    # Use the half distance in the first trace
    pix_01 = 1.5 * pols[0] - 0.5 * pols[1]

    borders.append((pix_01, pix_12))

    for idx in range(1, len(pols)-1):
        pix_01 = pix_12
        pix_12 = 0.5 * (pols[idx] + pols[idx+1])
        borders.append((pix_01, pix_12))

    # Estimate right border for the last trace
    pix_01 = pix_12
    # Use the half distance in last trace
    pix_12 = 1.5 * pols[-1] - 0.5 * pols[-2]

    borders.append((pix_01, pix_12))

    rss = extract_simple_rss(data, borders)

    return rss
