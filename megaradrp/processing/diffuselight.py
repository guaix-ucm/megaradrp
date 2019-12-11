#
# Copyright 2019 Universidad Complutense de Madrid
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


_logger = logging.getLogger(__name__)


class DiffuseLightCorrector(Corrector):
    """A Node that corrects a frame from diffuse light"""

    def __init__(self, diffuse, datamodel=None, calibid='calibid-unknown',
                 dtype='float32'):

        super(DiffuseLightCorrector, self).__init__(
            datamodel=datamodel,
            calibid=calibid,
            dtype=dtype)

        if isinstance(diffuse, fits.HDUList):
            self.corr = diffuse[0].data
        elif isinstance(diffuse, fits.ImageHDU):
            self.corr = diffuse.data
        else:
            self.corr = numpy.asarray(diffuse)

    def header_update(self, hdr, imgid):
        hdr['NUM-DFL'] = self.calibid
        hdr['history'] = 'Diffuse light correction {}'.format(imgid)
        hdr['history'] = 'Diffuse light correction time {}'.format(datetime.datetime.utcnow().isoformat())

    def run(self, img):
        imgid = self.get_imgid(img)
        _logger.debug('correct diffuse light in image %s', imgid)

        img['primary'].data -= self.corr
        hdr = img['primary'].header

        self.header_update(hdr, imgid)

        return img
