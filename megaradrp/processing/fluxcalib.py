#
# Copyright 2017 Universidad Complutense de Madrid
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
from numina.flow.processing import Corrector


_logger = logging.getLogger(__name__)


class FluxCalibration(Corrector):
    """A Node that calibrates absolute flux."""

    def __init__(self, sensitivity, datamodel=None, dtype='float32'):

        self.sensp = sensitivity['primary']
        self.limf = (self.sensp.header['PIXLIMF1'] - 1, self.sensp.header['PIXLIMF2'] - 1)
        self.limr = (self.sensp.header['PIXLIMR1'] - 1, self.sensp.header['PIXLIMR2'] - 1)
        self.limm = (self.sensp.header['PIXLIMM1'] - 1, self.sensp.header['PIXLIMM2'] - 1)

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
        newdata = numpy.zeros_like(img['primary'].data)
        validr = slice(self.limr[0], self.limr[1] + 1, 1)
        # with numpy.errstate(invalid='ignore', divide='ignore'):
        newdata[:, validr] = img['primary'].data[:, validr] / self.sensp.data[validr]
        img['primary'].data = newdata
        self.header_update(hdr, imgid)

        return img

    def header_update(self, hdr, imgid):
        hdr['BUNIT'] = self.target_unit
        hdr['NUM-STD'] = self.calibid
        hdr['NUM-FLUX'] = self.calibid
        hdr['PIXLIMF1'] = (self.limf[0] + 1, "Start of valid flux calibration")
        hdr['PIXLIMF2'] = (self.limf[1] + 1, "End of valid flux calibration")
        hdr['PIXLIMR1'] = self.limr[0] + 1
        hdr['PIXLIMR2'] = self.limr[1] + 1
        hdr['PIXLIMM1'] = self.limm[0] + 1
        hdr['PIXLIMM2'] = self.limm[1] + 1
        hdr['history'] = 'Flux calibration in {}'.format(imgid)
        hdr['history'] = 'Flux calibration time {}'.format(datetime.datetime.utcnow().isoformat())
