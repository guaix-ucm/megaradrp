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


def apextract_tracemap(data, tracemap):
    """Extract apertures using a tracemap.

    Consider that the nearest fiber could be far away if there
    are missing fibers.

    Parameters
    ----------

    data: ndarray
    tracemap: TraceMap

    """

    existing = [None]
    for t in tracemap.contents:
        if t.valid:
            existing.append(t)

    existing.append(None)
    # Compute borders
    borders = []
    far_dist = 100

    # Handle the first and last using centinels
    for t1, t2, t3 in zip(existing, existing[1:], existing[2:]):
        # Distance to contfibers # in box
        if t1 is None:
            d21 = far_dist
        else:
            d21 = (t2.fibid - t1.fibid) + (t2.boxid - t1.boxid)
        if t3 is None:
            d32 = far_dist
        else:
            d32 = (t3.fibid - t2.fibid) + (t3.boxid - t2.boxid)

        y2_off = t2.polynomial(tracemap.ref_column)
        t2_off = tracemap.global_offset(y2_off)
        t2_pol_off = t2.polynomial + t2_off

        # Right border
        pix_32 = None
        if d32 <= 2:
            y3_off = t3.polynomial(tracemap.ref_column)
            t3_off = tracemap.global_offset(y3_off)
            pix_32 = 0.5 * (t2_pol_off + t3.polynomial + t3_off)
            if d32 == 2:
                pix_32 = 0.5 * (pix_32 + t2_pol_off)
        elif d32 > 2:
            pix_32 = None

        pix_21 = None
        if d21 <= 2:
            y1_off = t1.polynomial(tracemap.ref_column)
            t1_off = tracemap.global_offset(y1_off)
            pix_21 = 0.5 * (t2_pol_off + t1.polynomial + t1_off)
            if d21 == 2:
                pix_21 = 0.5 * (pix_21 + t2_pol_off)
        elif d21 > 2:
            pix_21 = None

        if pix_32 is None and pix_21 is None:
            continue

        if pix_32 is None:
            # Recompute pix32 using pix_21
            pix_32 = t2_pol_off + (t2_pol_off - pix_21)

        if pix_21 is None:
            # Recompute pix21 using pix_32
            pix_21 = t2_pol_off - (pix_32 - t2_pol_off)

        borders.append((t2.fibid, pix_21, pix_32))

    nfibers = tracemap.total_fibers
    out = numpy.zeros((nfibers, data.shape[1]), dtype='float')
    rss = extract_simple_rss2(data, borders, out=out)

    return rss


def extract_simple_rss(arr, borders2, axis=0, out=None):

    # FIXME, this should be changed in numina
    # If arr is not in native byte order, the C-extension won't work
    if arr.dtype.byteorder != '=':
        arr2 = arr.byteswap().newbyteorder()
    else:
        arr2 = arr

    if axis == 0:
        arr3 = arr2
    elif axis == 1:
        arr3 = arr2.T
    else:
        raise ValueError("'axis' must be 0 or 1")

    if out is None:
        out = numpy.zeros((borders2[-1][0], arr3.shape[1]), dtype='float')

    xx = numpy.arange(arr3.shape[1])

    # Borders contains a list of function objects
    for idx, b1, b2 in borders2:
        bb1 = b1(xx)
        bb1[bb1 < -0.5] = -0.5
        bb2 = b2(xx)
        bb2[bb2 > arr3.shape[0] - 0.5] = arr3.shape[0] - 0.5
        extract.extract_simple_intl(arr3, xx, bb1, bb2, out[idx-1])
    return out


apextract_tracemap_2 = apextract_tracemap
extract_simple_rss2 = extract_simple_rss


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
