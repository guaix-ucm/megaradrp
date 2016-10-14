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

from __future__ import print_function

import numpy
import numpy.polynomial.polynomial as nppol

from numina.array.trace.extract import extract_simple_rss


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

    rss = extract_simple_rss(data, borders, out=out)

    return rss


def apextract_tracemap_2(data, tracemap):
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
    rss = extract_simple_rss(data, borders, out=out)

    return rss



