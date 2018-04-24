#
# Copyright 2011-2018 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import logging
import datetime

import numpy as np

from megaradrp.core.processing import trimOut, get_conf_value
from numina.processing import Corrector

_logger = logging.getLogger(__name__)


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

        fits.writeto('eq_estimado.fits', data, clobber=True)

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
