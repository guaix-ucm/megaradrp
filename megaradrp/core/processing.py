#
# Copyright 2011-2017 Universidad Complutense de Madrid
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

from __future__ import print_function, division

import numpy
import numpy.polynomial.polynomial as nppol
from astropy.io import fits

from numina.array.trace.extract import extract_simple_rss
from numina.array.trace.extract import extract_simple_intl

_direc = ['normal', 'mirror']


def trimOut(img, detconf, direction='normal', out='trimmed.fits'):
    '''

    :param data: original data to be trimmed. Can be either a string <path>, an ImageHDU or a numpy array
    :param direction: Can be either 'normal' or 'mirror'
    :param out: When the original data is a path, the returned value is a file too
    :param bins: is the index of the _binning dictionary
    :return:
    '''

    if issubclass(str, type(img)):
        with fits.open(img) as hdul:
            hdu = trim_and_o_array(hdul[0].data, detconf, direction=direction)
            fits.writeto(out, hdu, clobber=True)

    elif isinstance(img, fits.PrimaryHDU):
        finaldata = trim_and_o_array(img.data, detconf, direction=direction)
        img.data = finaldata
        return img

    elif isinstance(img, numpy.ndarray):
        return trim_and_o_array(img, detconf, direction=direction)


def trim_and_o_array(array, detconf, direction='normal'):
    """Trim a MEGARA array with overscan."""

    if direction not in _direc:
        raise ValueError("%s must be either 'normal' or 'mirror'" % direction)

    if direction == 'normal':
        direcfun = lambda x: x
    else:
        direcfun = numpy.fliplr

    trim1 = get_conf_value(detconf, 'trim1')
    trim2 = get_conf_value(detconf, 'trim2')
    bng = get_conf_value(detconf, 'bng')

    if bng not in [[1,1],[1,2],[2,1],[2,2]]:
        raise ValueError("%s must be one if '11', '12', '21, '22'" % bng)

    nr2 = ((trim1[0][1] - trim1[0][0]) + (trim2[0][1]-trim2[0][0])) // bng[0]
    nc2 = (trim1[1][1] - trim1[1][0]) // bng[1]
    nr = (trim1[0][1] - trim1[0][0]) // bng[0]

    trim1[0][0] //= bng[0]
    trim1[0][1] //= bng[0]
    trim2[0][0] //= bng[0]
    trim2[0][1] //= bng[0]

    trim1[1][0] //= bng[1]
    trim1[1][1] //= bng[1]
    trim2[1][0] //= bng[1]
    trim2[1][1] //= bng[1]

    finaldata = numpy.empty((nr2, nc2), dtype='float32')

    finY = trim1[0][1] + trim2[0][0]
    finaldata[:nr, :] = direcfun(array[:trim1[0][1], trim1[1][0]:trim1[1][1]])
    finaldata[nr:, :] = direcfun(array[trim2[0][0]:finY, trim2[1][0]:trim2[1][1]])

    return finaldata


def get_conf_value(confFile, key):
    return confFile[key]


def apextract(data, trace):
    """Extract apertures."""
    rss = numpy.empty((trace.shape[0], data.shape[1]), dtype='float32')
    for idx, r in enumerate(trace):
        l = r[0]
        r = r[2] + 1
        sl = (slice(l, r), )
        m = data[sl].sum(axis=0)
        rss[idx] = m
    return rss


def apextract_tracemap_2(data, tracemap):
    """Extract apertures using a tracemap.

    Consider that the nearest fiber could be far away if there
    are missing fibers.

    Parameters
    ----------

    data: ndarray
    tracemap: TraceMap

    """

    existing = [None]
    for t in tracemap.contents:
        if t.valid:
            existing.append(t)

    existing.append(None)
    # Compute borders
    borders = []
    far_dist = 100

    # Handle the first and last using centinels
    for t1, t2, t3 in zip(existing, existing[1:], existing[2:]):
        # Distance to contfibers # in box
        if t1 is None:
            d21 = far_dist
        else:
            d21 = (t2.fibid - t1.fibid) + (t2.boxid - t1.boxid)
        if t3 is None:
            d32 = far_dist
        else:
            d32 = (t3.fibid - t2.fibid) + (t3.boxid - t2.boxid)

        y2_off = t2.polynomial(tracemap.ref_column)
        t2_off = tracemap.global_offset(y2_off)
        t2_pol_off = t2.polynomial + t2_off

        # Right border
        pix_32 = None
        if d32 <= 2:
            y3_off = t3.polynomial(tracemap.ref_column)
            t3_off = tracemap.global_offset(y3_off)
            pix_32 = 0.5 * (t2_pol_off + t3.polynomial + t3_off)
            if d32 == 2:
                pix_32 = 0.5 * (pix_32 + t2_pol_off)
        elif d32 > 2:
            pix_32 = None

        pix_21 = None
        if d21 <= 2:
            y1_off = t1.polynomial(tracemap.ref_column)
            t1_off = tracemap.global_offset(y1_off)
            pix_21 = 0.5 * (t2_pol_off + t1.polynomial + t1_off)
            if d21 == 2:
                pix_21 = 0.5 * (pix_21 + t2_pol_off)
        elif d21 > 2:
            pix_21 = None

        if pix_32 is None and pix_21 is None:
            continue

        if pix_32 is None:
            # Recompute pix32 using pix_21
            pix_32 = t2_pol_off + (t2_pol_off - pix_21)

        if pix_21 is None:
            # Recompute pix21 using pix_32
            pix_21 = t2_pol_off - (pix_32 - t2_pol_off)

        borders.append((t2.fibid, pix_21, pix_32))

    nfibers = tracemap.total_fibers
    out = numpy.zeros((nfibers, data.shape[1]), dtype='float')
    rss = extract_simple_rss2(data, borders, out=out)

    return rss


def extract_simple_rss2(arr, borders2, axis=0, out=None):

    # FIXME, this should be changed in numina
    # If arr is not in native byte order, the C-extension won't work
    if arr.dtype.byteorder != '=':
        arr2 = arr.byteswap().newbyteorder()
    else:
        arr2 = arr

    if axis == 0:
        arr3 = arr2
    elif axis == 1:
        arr3 = arr2.T
    else:
        raise ValueError("'axis' must be 0 or 1")

    if out is None:
        out = numpy.zeros((borders2[-1][0], arr3.shape[1]), dtype='float')

    xx = numpy.arange(arr3.shape[1])

    # Borders contains a list of function objects
    for idx, b1, b2 in borders2:
        bb1 = b1(xx)
        bb1[bb1 < -0.5] = -0.5
        bb2 = b2(xx)
        bb2[bb2 > arr3.shape[0] - 0.5] = arr3.shape[0] - 0.5
        extract_simple_intl(arr3, xx, bb1, bb2, out[idx-1])
    return out


def apextract_tracemap(data, tracemap):
    """Extract apertures using a tracemap."""

    pols = [t.polynomial for t in tracemap.contents]

    borders = []

    # These are polynomial, they can be summed

    # Estimate left border for the first trace
    pix_12 = 0.5 * (pols[0] + pols[1])
    # Use the half distance in the first trace
    pix_01 = 1.5 * pols[0] - 0.5 * pols[1]

    borders.append((pix_01, pix_12))

    for idx in range(1, len(pols)-1):
        if pols[idx].order!=0:
            pix_01 = pix_12
            pix_12 = 0.5 * (pols[idx] + pols[idx+1])
            borders.append((pix_01, pix_12))
        else:
            empty = nppol.Polynomial([0.0])
            borders.append((empty, empty))
        # else:
        #     if pols[idx].order==0:
        #         borders.append((pix_01, pols[idx+1]))
        #     else:
        #         borders.append((pix_01, pols[idx+1]))

        # borders.append((pix_01, pix_12))
    # Estimate right border for the last trace
    pix_01 = pix_12
    # Use the half distance in last trace
    pix_12 = 1.5 * pols[-1] - 0.5 * pols[-2]

    borders.append((pix_01, pix_12))

    out = numpy.zeros((len(pols), data.shape[1]), dtype='float')

    rss = extract_simple_rss(data, borders, out=out)

    return rss


def apextract_weights(data, weights, size=4096):
    """Extract apertures using a map of weights."""

    import multiprocessing as mp

    def decompress(tar_file):
        '''
        :param tar_name: <str> name of the tar file
        :return: None
        '''

        name = tar_file.fileobj.name.split('.tar')[0]
        aux = tar_file.extractall(name+'/')
        return name


    processes = mp.cpu_count()-2

    path = decompress(weights)

    pool = mp.Pool(processes=processes)
    results = [pool.apply_async(load_files_paralell,
                                args=(ite, path)) for ite in range(size)]
    results = [p.get() for p in results]

    pool2 = mp.Pool(processes=processes)
    extracted_w = [pool2.apply_async(extract_w_paralell,
                                args=(data[:,ite], results[ite])) for ite in range(size)]
    extracted_w = [p.get() for p in extracted_w]

    return numpy.array(extracted_w)


def load_files_paralell(col, path):
    '''
    :param col: <str,int> name of the fits file. It is a counter
    :param path: <str> path where *.npz are
    :return: csr_matrix
    '''
    from scipy.sparse import csr_matrix

    filename = '%s/%s.npz' % (path, col)
    loader = numpy.load(filename)
    return csr_matrix(
        (loader['data'], loader['indices'], loader['indptr']),
        shape=loader['shape'])

from scipy.sparse.linalg import lsqr

def extract_w_paralell(img, mlist):
    '''
    :param img: <fits>
    :param mlist: <list> one element of the csr_matrix
    :return: <ndarray> result of lsqr
    '''
    x = lsqr(mlist, img)
    return x[0]
