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

from astropy.io import fits
import numpy
from numina.processing import Corrector
from numina.frame.utils import copy_img


_logger = logging.getLogger(__name__)


class CommonFlatCorrector(Corrector):
    """A Node that corrects from fiber flat."""

    def __init__(self, flatdata, datamodel=None, calibid='calibid-unknown', dtype='float32'):

        super(CommonFlatCorrector, self).__init__(datamodel=datamodel,
                                                 calibid=calibid,
                                                 dtype=dtype)
        if isinstance(flatdata, fits.HDUList):
            self.corr = flatdata[0].data
        elif isinstance(flatdata, fits.ImageHDU):
            self.corr = flatdata.data
        else:
            self.corr = numpy.asarray(flatdata)

        self.corrmean = self.corr.mean()
        self.flattag = 'flat'

    def run(self, img):
        imgid = self.get_imgid(img)
        _logger.debug('correct %s in image %s', self.flattag, imgid)

        # Avoid nan values when divide
        my_mask = self.corr == 0.0
        self.corr[my_mask] = 1.0

        img['primary'].data /= self.corr
        hdr = img['primary'].header

        self.header_update(hdr, imgid)

        return img

    def header_update(self, hdr, imgid):
        hdr['NUM-FLT'] = self.calibid
        hdr['history'] = 'Flat correction {}'.format(imgid)
        hdr['history'] = 'Flat correction time {}'.format(datetime.datetime.utcnow().isoformat())


class FiberFlatCorrector(CommonFlatCorrector):
    """A Node that corrects from fiber flat."""

    def __init__(self, fiberflat, datamodel=None, dtype='float32'):
        calibid = datamodel.get_imgid(fiberflat)
        super(FiberFlatCorrector, self).__init__(
            flatdata=fiberflat,
            datamodel=datamodel,
            calibid=calibid,
            dtype=dtype
        )
        self.flattag = 'fiberflat'

    def header_update(self, hdr, imgid):
        hdr['NUM-FIBF'] = self.calibid
        hdr['history'] = 'Fiber flat correction {}'.format(imgid)
        hdr['history'] = 'Fiber flat correction time {}'.format(datetime.datetime.utcnow().isoformat())


class Splitter(Corrector):
    """Make a copy of its input"""
    def __init__(self):
        super(Splitter, self).__init__()
        self.out = None

    def run(self, img):
        self.out = self.copy_img(img)
        return img

    def copy_img(self, img):
        return copy_img(img)


class FlipLR(Corrector):
    """Make a copy of its input"""
    def __init__(self):
        super(FlipLR, self).__init__()

    def run(self, img):
        img[0].data = numpy.fliplr(img[0].data)
        return img
