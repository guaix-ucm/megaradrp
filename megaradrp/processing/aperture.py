#
# Copyright 2011-2017 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import logging
import datetime

import astropy.io.fits as fits
from numina.flow.processing import Corrector

from megaradrp.core.processing import apextract_tracemap_2

_logger = logging.getLogger(__name__)


class ApertureExtractor(Corrector):
    """A Node that extracts apertures."""

    def __init__(self, tracemap, datamodel=None, dtype='float32', processes=0):

        self.tracemap = tracemap
        self.processes = processes
        super(ApertureExtractor, self).__init__(datamodel=datamodel,
                                                calibid=tracemap.uuid,
                                                dtype=dtype)

    def run(self, img):
        # workaround
        imgid = self.get_imgid(img)

        method_name = 'simple'
        simple = True
        if hasattr(self.tracemap, 'aper_extract'):
            simple = False
            method_name = 'advanced'

        if simple:
            _logger.debug('simple aperture extraction')
            _logger.debug('extracting (apextract_tracemap) in image %s', imgid)
            _logger.debug('with trace map %s', self.calibid)
        else:
            _logger.debug('advanced aperture extraction')
            _logger.debug('extracting (apextract_model) in image %s', imgid)
            _logger.debug('with model map %s', self.calibid)

        _logger.debug('offsets are %s', self.tracemap.global_offset.coef)
        if simple:
            rssdata = apextract_tracemap_2(img[0].data, self.tracemap)
        else:
            rssdata = self.tracemap.aper_extract(img[0].data, processes=self.processes)

        img[0].data = rssdata

        hdr = img[0].header

        hdr['NUM-APE'] = self.calibid
        hdr['history'] = 'Aperture extraction method {}'.format(method_name)
        hdr['history'] = 'Aperture extraction with {}'.format(self.calibid)
        hdr['history'] = 'Aperture extraction offsets are {}'.format(self.tracemap.global_offset.coef.tolist())
        hdr['history'] = 'Aperture extraction time {}'.format(datetime.datetime.utcnow().isoformat())

        # Update Fibers
        fibers_ext = img['FIBERS']
        fibers_ext_headers = fibers_ext.header
        for aper in self.tracemap.contents:
            key = "FIB%03d_V" % aper.fibid
            fibers_ext_headers[key] = aper.valid
            key = "FIB%03dS1" % aper.fibid
            fibers_ext_headers[key] = aper.start
            key = "FIB%03dS2" % aper.fibid
            fibers_ext_headers[key] = aper.stop

        newimg = fits.HDUList([img[0], fibers_ext])
        return newimg

