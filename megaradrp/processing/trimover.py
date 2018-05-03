#
# Copyright 2011-2018 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

from __future__ import print_function, division

import logging
import datetime

import numpy as np
import numpy
import numpy.polynomial.polynomial as nppol
from astropy.io import fits
import numina.array.trace.extract as extract
from numina.processing import Corrector


_logger = logging.getLogger(__name__)


_direc = ['normal', 'mirror']


def trimOut(img, detconf, direction='normal', out='trimmed.fits'):
    """
    :param data: original data to be trimmed. Can be either a string <path>, an ImageHDU or a numpy array
    :param direction: Can be either 'normal' or 'mirror'
    :param out: When the original data is a path, the returned value is a file too
    :param bins: is the index of the _binning dictionary
    :return:
    """

    if issubclass(str, type(img)):
        with fits.open(img) as hdul:
            hdu = trim_and_o_array(hdul[0].data, detconf, direction=direction)
            fits.writeto(out, hdu, overwrite=True)

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


def apextract_weights(data, weights, size=4096):
    """Extract apertures using a map of weights."""

    import multiprocessing as mp

    def decompress(tar_file):
        """
        :param tar_name: <str> name of the tar file
        :return: None
        """

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
    """"
    :param col: <str,int> name of the fits file. It is a counter
    :param path: <str> path where *.npz are
    :return: csr_matrix
    """
    from scipy.sparse import csr_matrix

    filename = '%s/%s.npz' % (path, col)
    loader = numpy.load(filename)
    return csr_matrix(
        (loader['data'], loader['indices'], loader['indptr']),
        shape=loader['shape'])



def extract_w_paralell(img, mlist):
    """
    :param img: <fits>
    :param mlist: <list> one element of the csr_matrix
    :return: <ndarray> result of lsqr
    """

    from scipy.sparse.linalg import lsqr
    x = lsqr(mlist, img)
    return x[0]


class OverscanCorrector(Corrector):
    """A Corrector Node that corrects a MEGARA image from overscan."""

    def __init__(self, detconf, datamodel=None, calibid='calibid-unknown', dtype='float32'):

        trim1 = get_conf_value(detconf, 'trim1')
        trim2 = get_conf_value(detconf, 'trim2')
        bng = get_conf_value(detconf, 'bng')
        overscan1 = get_conf_value(detconf, 'overscan1')
        overscan2 = get_conf_value(detconf, 'overscan2')
        prescan1 = get_conf_value(detconf, 'prescan1')
        prescan2 = get_conf_value(detconf, 'prescan2')
        middle1 = get_conf_value(detconf, 'middle1')
        middle2 = get_conf_value(detconf, 'middle2')

        auxX, auxY, auxZ, auxT = self.data_binning(trim1, bng)
        middleX, middleY, middleZ, middleT = self.data_binning(middle1, bng)
        prescanX, prescanY, prescanZ, prescanT = self.data_binning(prescan1,
                                                                   bng)
        overscanX, overscanY, overscanZ, overscanT = self.data_binning(
            overscan1, bng)
        self.trim1 = (slice(auxX, auxY), slice(auxZ, auxT))
        self.orow1 = (slice(middleX, middleY), slice(middleZ, middleT))
        self.pcol1 = (slice(prescanX, prescanY), slice(prescanZ, prescanT))
        self.ocol1 = (slice(overscanX, overscanY), slice(overscanZ, overscanT))

        auxX, auxY, auxZ, auxT = self.data_binning(trim2, bng)
        middleX, middleY, middleZ, middleT = self.data_binning(middle2, bng)
        prescanX, prescanY, prescanZ, prescanT = self.data_binning(prescan2,
                                                                   bng)
        overscanX, overscanY, overscanZ, overscanT = self.data_binning(
            overscan2, bng)
        self.trim2 = (slice(auxX, auxY), slice(auxZ, auxT))
        self.orow2 = (slice(middleX, middleY), slice(middleZ, middleT))
        self.pcol2 = (slice(prescanX, prescanY), slice(prescanZ, prescanT))
        self.ocol2 = (slice(overscanX, overscanY), slice(overscanZ, overscanT))

        # self.test_image()
        super(OverscanCorrector, self).__init__(datamodel=datamodel,
                                                calibid=calibid,
                                                dtype=dtype)

    def data_binning(self, data, binning):
        '''
         Axis:x --> factorX
         Axis:y --> factorY
        '''
        factorX = 1.0 / binning[1]
        factorY = 1.0 / binning[0]
        x = int(factorY * data[0][0])
        y = int(factorY * data[0][1])
        z = int(factorX * data[1][0])
        t = int(factorX * data[1][1])

        return x, y, z, t

    def test_image(self):
        import astropy.io.fits as fits

        data = np.empty((4212, 4196), dtype='float32')
        data[self.pcol1] += 1
        data[self.orow1] += 10
        data[self.ocol1] += 100
        data[self.trim1] += 1000

        data[self.pcol2] += 5
        data[self.orow2] += 50
        data[self.ocol2] += 500
        data[self.trim2] += 5000

        fits.writeto('eq_estimado.fits', data, overwrite=True)

    def run(self, img):
        imgid = self.get_imgid(img)
        data = img[0].data

        p1 = data[self.pcol1].mean()
        _logger.debug('prescan1 is %f', p1)
        # or1 = data[self.orow1].mean()
        # _logger.debug('row overscan1 is %f', or1)
        oc1 = data[self.ocol1].mean()
        _logger.debug('col overscan1 is %f', oc1)
        # avg = (p1 + or1 + oc1) / 3.0
        avg = (p1 + oc1) / 2.0
        _logger.debug('average scan1 is %f', avg)
        data[self.trim1] -= avg

        p2 = data[self.pcol2].mean()
        _logger.debug('prescan2 is %f', p2)
        # or2 = data[self.orow2].mean()
        # _logger.debug('row overscan2 is %f', or2)
        oc2 = data[self.ocol2].mean()
        _logger.debug('col overscan2 is %f', oc2)
        # avg = (p2 + or2 + oc2) / 3.0
        avg = (p2 + oc2) / 2.0
        _logger.debug('average scan2 is %f', avg)
        data[self.trim2] -= avg
        hdr = img['primary'].header
        hdr['NUM-OVPE'] = self.calibid
        hdr['history'] = 'Overscan correction with {}'.format(self.calibid)
        hdr['history'] = 'Overscan correction time {}'.format(datetime.datetime.utcnow().isoformat())
        hdr['history'] = 'Mean of prescan1 is %f' % p1
        hdr['history'] = 'col overscan1 is %f' %  oc1
        hdr['history'] = 'average scan1 is %f' % avg
        hdr['history'] = 'prescan2 is %f' %  p2
        hdr['history'] = 'col overscan2 is %f' % oc2
        hdr['history'] = 'average scan2 is %f' % avg

        return img


class TrimImage(Corrector):
    """A Corrector Node that trims MEGARA images."""

    def __init__(self, detconf, datamodel=None, calibid='calibid-unknown', dtype='float32'):
        self.detconf = detconf
        super(TrimImage, self).__init__(
            datamodel=datamodel,
            calibid=calibid,
            dtype=dtype
        )

    def run(self, img):
        imgid = self.get_imgid(img)
        _logger.debug('trimming image %s', imgid)

        img[0] = trimOut(img[0], self.detconf)
        hdr = img['primary'].header
        hdr['NUM-TRIM'] = self.calibid
        hdr['history'] = 'Trimming correction with {}'.format(self.calibid)
        hdr['history'] = 'Trimming correction time {}'.format(datetime.datetime.utcnow().isoformat())

        return img


class GainCorrector(Corrector):
    """A Corrector Node that corrects MEGARA images from different gain."""

    def __init__(self, detconf, datamodel=None, calibid='calibid-unknown', dtype='float32'):
        self.detconf = detconf
        self.gain1 = self.detconf['gain1']
        self.gain2 = self.detconf['gain2']
        super(GainCorrector, self).__init__(
            datamodel=datamodel,
            calibid=calibid,
            dtype=dtype
        )

    def run(self, img):
        imgid = self.get_imgid(img)
        _logger.debug('gain correction in image %s', imgid)

        # img[0] = trimOut(img[0], self.detconf)
        hdr = img['primary'].header
        hdr['NUM-GAIN'] = self.calibid
        hdr['BUNIT'] = 'ELECTRON'

        part = img[0].data.shape[0] // 2

        img[0].data[:part] *= self.gain1
        img[0].data[part:] *= self.gain2

        hdr['history'] = 'Gain correction with {}'.format(self.calibid)
        hdr['history'] = 'Gain correction time {}'.format(datetime.datetime.utcnow().isoformat())
        hdr['history'] = 'Gain1 correction value {}'.format(self.gain1)
        hdr['history'] = 'Gain2 correction value {}'.format(self.gain2)



        return img
