#
# Copyright 2017-2019 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import logging
import datetime


import astropy.wcs
import astropy.units as u

import numpy
from numina.processing import Corrector


_logger = logging.getLogger(__name__)


PIXLIM_KEYS = ['PIXLIMF1', 'PIXLIMF2', 'PIXLIMR1', 'PIXLIMR2', 'PIXLIMM1', 'PIXLIMM2']
WAVLIM_KEYS = ['WAVLIMF1', 'WAVLIMF2', 'WAVLIMR1', 'WAVLIMR2', 'WAVLIMM1', 'WAVLIMM2']


class FluxCalibration(Corrector):
    """A Node that calibrates absolute flux."""

    def __init__(self, sensitivity, datamodel=None, dtype='float32'):

        self.sensp = sensitivity['primary']
        self.pixlims = {key: (self.sensp.header[key] - 1) for key in PIXLIM_KEYS}

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
        limr1 = self.pixlims['PIXLIMR1']
        limr2 = self.pixlims['PIXLIMR2']
        validr = slice(limr1, limr2 + 1, 1)
        # with numpy.errstate(invalid='ignore', divide='ignore'):
        newdata[:, validr] = img['primary'].data[:, validr] / self.sensp.data[validr]
        img['primary'].data = newdata
        self.header_update(hdr, imgid)

        return img

    def header_update(self, hdr, imgid):
        hdr['BUNIT'] = self.target_unit
        hdr['NUM-STD'] = self.calibid
        hdr['NUM-FLUX'] = self.calibid

        update_flux_limits(hdr, self.pixlims, ref=0)

        hdr['history'] = 'Flux calibration in {}'.format(imgid)
        hdr['history'] = 'Flux calibration time {}'.format(datetime.datetime.utcnow().isoformat())


def update_flux_limits(header, pixlims, wcs=None, ref=1):
    """Update keywords used for flux limits"""

    pixlim_coms = ["Start of valid flux calibration",
                   "End of valid flux calibration",
                   'Start of region with at least one fiber',
                   'End of region with at least one fiber',
                   'Start of region with all fibers',
                   'End of region with all fibers',
                   ]

    if ref not in [0, 1]:
        raise ValueError("ref must be 0 or 1")

    off = (ref + 1) % 2

    if wcs is None:
        wcs = astropy.wcs.WCS(header)

    for key, com in zip(PIXLIM_KEYS, pixlim_coms):
        header[key] = (pixlims[key] + off, com)

    r1 = numpy.array([pixlims[key] for key in PIXLIM_KEYS])
    r2 = numpy.zeros_like(r1)
    lm = numpy.array([r1, r2])
    # Values are 0-based
    wavelen_ = wcs.all_pix2world(lm.T, ref)
    if wcs.wcs.cunit[0] == u.dimensionless_unscaled:
        # CUNIT is empty, assume Angstroms
        wavelen = wavelen_[:, 0] * u.AA
    else:
        wavelen = wavelen_[:, 0] * wcs.wcs.cunit[0]

    for idx, (key, com) in enumerate(zip(WAVLIM_KEYS, pixlim_coms)):
        header[key] = (wavelen[idx].to(u.AA).value, com)

    return header