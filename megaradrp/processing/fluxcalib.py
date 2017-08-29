#
# Copyright 2017 Universidad Complutense de Madrid
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


_logger = logging.getLogger(__name__)


class FluxCalibration(Corrector):
    """A Node that calibrates absolute flux."""

    def __init__(self, sensitivity, datamodel=None, dtype='float32'):

        self.sensp = sensitivity['primary']
        calibid = datamodel.get_imgid(sensitivity)
        self.target_unit = self.sensp.header['TUNIT']
        super(FluxCalibration, self).__init__(datamodel=datamodel,
                                              calibid=calibid,
                                              dtype=dtype)

    def run(self, img):
        imgid = self.get_imgid(img)
        _logger.debug('calibrate flux in image %s', imgid)
        hdr = img['primary'].header
        exptime = hdr['EXPTIME']
        img['primary'].data /= exptime
        img['primary'].data *= self.sensp.data
        self.header_update(hdr, imgid)

        return img

    def header_update(self, hdr, imgid):
        hdr['BUNIT'] = self.target_unit
        hdr['NUM-STD'] = self.calibid
        hdr['NUM-FLUX'] = self.calibid
        hdr['history'] = 'Flux calibration in {}'.format(imgid)
        hdr['history'] = 'Flux calibration time {}'.format(datetime.datetime.utcnow().isoformat())
