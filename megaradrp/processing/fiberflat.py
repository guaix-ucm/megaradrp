#
# Copyright 2011-2016 Universidad Complutense de Madrid
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

import logging
import datetime

from astropy.io import fits
import numpy
from numina.flow.processing import Corrector

_logger = logging.getLogger('numina.processing')


class FiberFlatCorrector(Corrector):
    '''A Node that corrects from fiber flat.'''

    def __init__(self, fiberflat, datamodel=None, dtype='float32'):

        super(FiberFlatCorrector, self).__init__(datamodel=datamodel,
                                                 dtype=dtype)
        if isinstance(fiberflat, fits.HDUList):
            self.corr = fiberflat[0].data
        elif isinstance(fiberflat, fits.ImageHDU):
            self.corr = fiberflat.data
        else:
            self.corr = numpy.asarray(fiberflat)

        self.corrmean = self.corr.mean()
        # self.corrid = self.get_imgid(fiberflat)

    def run(self, img):
        imgid = self.get_imgid(img)
        _logger.debug('correct from fiber flat in image %s', imgid)

        # Avoid nan values when divide
        my_mask = self.corr == 0.0
        self.corr[my_mask] = 1.0

        img['primary'].data /= self.corr
        hdr = img['primary'].header

        self.header_update(hdr, imgid)

        return img

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
        return fits.HDUList([hdu.copy() for hdu in img])


class FlipLR(Corrector):
    """Make a copy of its input"""
    def __init__(self):
        super(FlipLR, self).__init__()

    def run(self, img):
        img[0].data = numpy.fliplr(img[0].data)
        return img
