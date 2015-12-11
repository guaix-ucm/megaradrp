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
        hdu = trim_and_o_hdu(hdul[0])
        hdu.writeto(out, clobber=True)


def trim_and_o_hdu(hdu):
    """Trim a MEGARA HDU with overscan."""

    # FIXME: this should come from the header
    direction = 'normal'
    bins = '11'

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

    # finaldata = np.empty((shape_number, shape_number), dtype='float32')
    # finaldata[0:shape_number//2,0:shape_number] = array[0:shape_number//2,PSCANW:shape_number+PSCANW]
    # finaldata[shape_number//2:shape_number,0:shape_number] = array[shape_number//2+2*OSCANW:shape_number+2*OSCANW,PSCANW:shape_number+PSCANW]

    finaldata = np.empty((nr2, nc2), dtype='float32')
    finaldata[:nr, :] = direcfun(array[:nr, psc1:nc2 + psc1])
    # finaldata[:nr, :] = direcfun(array[:nr,psc1:nc2+psc1])
    finaldata[nr:, :] = direcfun(array[nr + oscan2:, psc1:nc2 + psc1])
    # finaldata[nr:, :] = direcfun(array[nr+2*oscan2:nc2+2*oscan2, psc1:nc2 + psc1])
    # finaldata[nr:, :] = direcfun(array[nr+oscan2:nc2+2*oscan2, psc1:nc2 + psc1])
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


def fill_other(data, a, b):
    start = wcs_to_pix_1d(a)
    end = wcs_to_pix_1d(b)
    data[start] = min(start+0.5, b)-a
    data[start+1:end] = 1.0
    if end > start:
        data[end] = b - (end-0.5)
    return data


def extract_region(data, border1, border2, pesos, xpos):
        
    extend = (border1.min(), border2.max())
    extend_pix = (wcs_to_pix_1d(extend[0]), wcs_to_pix_1d(extend[1])+1)
    region = slice(extend_pix[0],extend_pix[1])

    for x, a,b in zip(xpos, border1, border2):
        fill_other(pesos[:,x], a, b)

    final2d = data[region,:] * pesos[region,:]

    pesos[region,:] = 0.0
    final = final2d.sum(axis=0)
    return final


def apextract_tracemap(data, tracemap):
    """Extract apertures using a tracemap."""
    
    # FIXME: a little hackish
    
    pols = [np.poly1d(t['fitparms']) for t in tracemap]
    
    borders = []
    pix_12 = 0.5 * (pols[1] + pols[0])
    # Use the half distance in the first trace
    pix_01 = 1.5 * pols[0] - 0.5 * pols[1]
    pix_2 = pols[1]    
    borders.append((pix_01, pix_12))

    for p2 in pols[2:-1]:
     
        pix_1, pix_01 = pix_2, pix_12 
        pix_2 = p2

        pix_12 = 0.5 * (pix_2 + pix_1)
        borders.append((pix_01, pix_12))

    # extract the last trace
    pix_1, pix_01 = pix_2, pix_12 
    pix_12 = 2 * pix_1 - pix_01
    borders.append((pix_01, pix_12))

    rss = extract_simple_rss(data, borders)

    return rss


