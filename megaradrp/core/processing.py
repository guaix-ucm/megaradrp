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
import numpy as np

from astropy.io import fits
from numina.array.trace.extract import extract_simple_rss

_direc = ['normal', 'mirror']

def trimOut(data, direction='normal', out='trimmed.fits', confFile={}):
    '''

    :param data: original data to be trimmed. Can be either a string <path>, an ImageHDU or a numpy array
    :param direction: Can be either 'normal' or 'mirror'
    :param out: When the original data is a path, the returned value is a file too
    :param bins: is the index of the _binning dictionary
    :return:
    '''

    if issubclass(str,type(data)):
        with fits.open(data) as hdul:
            hdu = trim_and_o_array(hdul[0].data, direction=direction, confFile=confFile)
            fits.writeto(out, hdu, clobber=True)

    elif isinstance(data,fits.PrimaryHDU):
        finaldata = trim_and_o_array(data.data, direction=direction, confFile=confFile)
        data.data = finaldata
        return data

    elif isinstance(data,np.ndarray):
        return trim_and_o_array(data, direction=direction, confFile=confFile)


def trim_and_o_array(array, direction='normal', confFile={}):
    """Trim a MEGARA array with overscan."""

    if direction not in _direc:
        raise ValueError("%s must be either 'normal' or 'mirror'" % direction)

    if direction == 'normal':
        direcfun = lambda x: x
    else:
        direcfun = np.fliplr

    trim1 = get_conf_value(confFile, 'trim1')
    trim2 = get_conf_value(confFile, 'trim2')
    bng = get_conf_value(confFile, 'bng')

    if bng not in [[1,1],[1,2],[2,1],[2,2]]:
        raise ValueError("%s must be one if '11', '12', '21, '22'" % bng)

    nr2 = ((trim1[0][1] - trim1[0][0]) + (trim2[0][1]-trim2[0][0]))/bng[0]
    nc2 = (trim1[1][1] - trim1[1][0]) / bng[1]
    nr = (trim1[0][1] - trim1[0][0])/bng[0]

    trim1[0][0] /= bng[0]
    trim1[0][1] /= bng[0]
    trim2[0][0] /= bng[0]
    trim2[0][1] /= bng[0]

    trim1[1][0] /= bng[1]
    trim1[1][1] /= bng[1]
    trim2[1][0] /= bng[1]
    trim2[1][1] /= bng[1]

    finaldata = np.empty((nr2, nc2), dtype='float32')

    finY = trim1[0][1] + trim2[0][0]
    finaldata[:nr, :] = direcfun(array[:trim1[0][1], trim1[1][0]:trim1[1][1]])
    finaldata[nr:, :] = direcfun(array[trim2[0][0]:finY, trim2[1][0]:trim2[1][1]])

    return finaldata

def get_conf_value(confFile, key=''):
    if confFile:
        if key in confFile.keys():
            if confFile[key]:
                return confFile[key]
            else:
                raise ValueError('Value not defined')
        else:
            raise ValueError('Key is not in configuration file')
    raise ValueError('Instrument configuration is not in the system')

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





def apextract_weights(data, weights, size=4096):
    """Extract apertures using a tracemap."""

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

    return np.array(extracted_w)


def load_files_paralell(col, path):
    '''
    :param col: <str,int> name of the fits file. It is a counter
    :param path: <str> path where *.npz are
    :return: csr_matrix
    '''
    from scipy.sparse import csr_matrix

    filename = '%s/%s.npz' % (path, col)
    loader = np.load(filename)
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