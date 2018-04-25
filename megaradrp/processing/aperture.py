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

import numpy
import numpy.polynomial.polynomial as nppol
import astropy.io.fits as fits
import numina.array.trace.extract as extract
import numina.processing


_logger = logging.getLogger(__name__)


def apextract(data, trace):
    """Extract apertures."""
    rss = numpy.empty((trace.shape[0], data.shape[1]), dtype='float32')
    for idx, r in enumerate(trace):
        l = r[0]
        r = r[2] + 1
        sl = (slice(l, r), )
        m = data[sl].sum(axis=0)
        rss[idx] = m
    return rss


def apextract_tracemap_old(data, tracemap):
    """Extract apertures using a tracemap."""

    pols = [t.polynomial for t in tracemap.contents]

    borders = []

    # These are polynomial, they can be summed

    # Estimate left border for the first trace
    pix_12 = 0.5 * (pols[0] + pols[1])
    # Use the half distance in the first trace
    pix_01 = 1.5 * pols[0] - 0.5 * pols[1]

    borders.append((pix_01, pix_12))

    for idx in range(1, len(pols)-1):
        if pols[idx].order!=0:
            pix_01 = pix_12
            pix_12 = 0.5 * (pols[idx] + pols[idx+1])
            borders.append((pix_01, pix_12))
        else:
            empty = nppol.Polynomial([0.0])
            borders.append((empty, empty))
        # else:
        #     if pols[idx].order==0:
        #         borders.append((pix_01, pols[idx+1]))
        #     else:
        #         borders.append((pix_01, pols[idx+1]))

        # borders.append((pix_01, pix_12))
    # Estimate right border for the last trace
    pix_01 = pix_12
    # Use the half distance in last trace
    pix_12 = 1.5 * pols[-1] - 0.5 * pols[-2]

    borders.append((pix_01, pix_12))

    out = numpy.zeros((len(pols), data.shape[1]), dtype='float')

    rss = extract.extract_simple_rss(data, borders, out=out)

    return rss


def apextract_tracemap(data, tracemap):
    """Extract apertures using a tracemap.

    Consider that the nearest fiber could be far away if there
    are missing fibers.
    """

    existing = [None]
    for t in tracemap.contents:
        if t.fitparms: # This is a check if the trace is valid
            existing.append(t)

    existing.append(None)
    # Compute borders
    borders = []

    # Handle the first and last using centinels
    for t1, t2, t3 in zip(existing, existing[1:], existing[2:]):
        # Distance to contfibers # in box
        if t1 is None:
            d21 = 100
        else:
            d21 = (t2.fibid - t1.fibid) + (t2.boxid - t1.boxid)
        if t3 is None:
            d32 = 100
        else:
            d32 = (t3.fibid - t2.fibid) + (t3.boxid - t2.boxid)

        # Right border
        if d32 == 1:
            pix_32 = 0.5 * (t2.polynomial + t3.polynomial)
        elif d32 == 2:
            pix_32 = 0.5 * (t2.polynomial + t3.polynomial)
            pix_32 = 0.5 * (pix_32 + t2.polynomial)
        elif d32 > 2:
            pix_32 = None

        if d21 == 1:
            pix_21 = 0.5 * (t2.polynomial + t1.polynomial)
        elif d21 == 2:
            pix_21 = 0.5 * (t2.polynomial + t1.polynomial)
            pix_21 = 0.5 * (pix_21 + t2.polynomial)
        elif d21 > 2:
            pix_21 = None

        if pix_32 is None and pix_21 is None:
            continue

        if pix_32 is None:
            # Recompute pix32 using pix_21
            pix_32 = t2.polynomial + (t2.polynomial - pix_21)

        if pix_21 is None:
            # Recompute pix21 using pix_32
            pix_21 = t2.polynomial - (pix_32 - t2.polynomial)

        #
        borders.append((t2.fibid - 1, (pix_21, pix_32)))

    mm = 623 # FIXME, hardcoded
    out = numpy.zeros((mm, data.shape[1]), dtype='float32')
    rss = extract.extract_simple_rss(data, borders, out=out)

    return rss


class ApertureExtractor(numina.processing.Corrector):
    """A Node that extracts apertures."""

    def __init__(self, trace_repr, datamodel=None, dtype='float32',
                 processes=0, offset=None):

        if offset:
            trace_repr.global_offset = trace_repr.global_offset + nppol.Polynomial(offset)

        self.trace_repr = trace_repr
        self.processes = processes
        super(ApertureExtractor, self).__init__(
            datamodel=datamodel,
            calibid=trace_repr.uuid,
            dtype=dtype
        )

    def run(self, img):
        # workaround
        imgid = self.get_imgid(img)

        method_name = 'simple'
        simple = True
        if hasattr(self.trace_repr, 'aper_extract'):
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

        _logger.debug('offsets are %s', self.trace_repr.global_offset.coef)
        if simple:
            rssdata = apextract_tracemap(img[0].data, self.trace_repr)
        else:
            rssdata = self.trace_repr.aper_extract(img[0].data, processes=self.processes)

        img[0].data = rssdata

        hdr = img[0].header

        hdr['NUM-APE'] = self.calibid
        hdr['history'] = 'Aperture extraction method {}'.format(method_name)
        hdr['history'] = 'Aperture extraction with {}'.format(self.calibid)
        hdr['history'] = 'Aperture extraction offsets are {}'.format(self.trace_repr.global_offset.coef.tolist())
        hdr['history'] = 'Aperture extraction time {}'.format(datetime.datetime.utcnow().isoformat())

        # Update Fibers
        fibers_ext = img['FIBERS']
        fibers_ext_headers = fibers_ext.header
        for aper in self.trace_repr.contents:
            key = "FIB%03d_V" % aper.fibid
            fibers_ext_headers[key] = aper.valid
            key = "FIB%03dS1" % aper.fibid
            fibers_ext_headers[key] = aper.start
            key = "FIB%03dS2" % aper.fibid
            fibers_ext_headers[key] = aper.stop

        newimg = fits.HDUList([img[0], fibers_ext])
        return newimg
