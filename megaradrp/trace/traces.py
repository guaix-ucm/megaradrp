#
# Copyright 2015 Universidad Complutense de Madrid
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

from __future__ import division, print_function

import numpy as np
import logging

from .peakdetection import peak_detection_mean_window  #Should be avoided
from numina.array.peaks.peakdet import refine_peaks
from skimage.feature import peak_local_max


_logger = logging.getLogger('numina.recipes.megara')

# Number of fibers in the boxes of the pseudo-slit of the LCB
boxes = [
    {'nfibers': 21,
     'id': 1},
    {'nfibers': 21,
     'id': 2},
    {'nfibers': 21,
     'id': 3},
    {'nfibers': 28,
     'id': 4},
    {'nfibers': 28,
     'id': 5},
    {'nfibers': 35,
     'id': 6},
    {'nfibers': 42,
     'id': 7},
    {'nfibers': 77,
     'id': 8},
    {'nfibers': 77,
     'id': 9},
    {'nfibers': 77,
     'id': 10},
    {'nfibers': 42,
     'id': 11},
    {'nfibers': 35,
     'id': 12},
    {'nfibers': 28,
     'id': 13},
    {'nfibers': 28,
     'id': 14},
    {'nfibers': 21,
     'id': 15},
    {'nfibers': 21,
     'id': 16},
    {'nfibers': 21,
     'id': 17}
]


def delicate_centre(x, y):
    pos = np.polyfit(x, y, deg=2)
    tx = -pos[1] / (2 * pos[0])
    py = np.polyval(pos, tx)
    return tx, py, pos


class FiberTraceInfo(object):
    def __init__(self, fibid, boxid):
        self.boxid = boxid
        self.fibid = fibid
        self.start = None


def init_traces(image, center, hs, background, maxdis=9.0):
    ixmin = 0
    ixmax = image.shape[0]

    xx = np.arange(image.shape[0])

    cut_region = slice(center - hs, center + hs)
    cut = image[:, cut_region]
    colcut = cut.mean(axis=1)
    maxt = peak_detection_mean_window(colcut, x=xx, k=3, xmin=ixmin, xmax=ixmax, background=background)

    peakdist = np.diff(maxt[:, 1])
    npeaks = maxt.shape[0]

    # Count the number of groups
    # fibers separated by more than "maxdis" pixels
    group_changes = np.nonzero(peakdist > maxdis)[0]
    # Fill groups
    groups = np.ones(npeaks, dtype='int')
    for p in group_changes:
        groups[p + 1:] += 1

    fiber_traces = {}
    for fibid in range(1, npeaks + 1):
        fiber_traces[fibid] = FiberTraceInfo(fibid, int(groups[fibid - 1]))

    # Window to fix the center of the trace
    fw = 2
    for fibid, trace in fiber_traces.items():
        _xi, _yi, _vali = center, maxt[fibid - 1, 1], maxt[fibid - 1, 2]
        pixmax = int(maxt[fibid - 1, 0])
        # Take 2*2+1 pix
        # This part and interp_max_3(image[nearp3-1:nearp3+2, col])
        # should do the same
        tx, py, _pos = delicate_centre(xx[pixmax - fw: pixmax + fw + 1],
                                       colcut[pixmax - fw: pixmax + fw + 1])

        trace.start = (center, tx, py)

    return fiber_traces


def init_traces_ex(image, center, hs, box_borders, tol=1.5):

    cut_region = slice(center-hs, center+hs)
    cut = image[:,cut_region]
    colcut = cut.mean(axis=1)

    counted_fibers = 0
    fiber_traces = {}
    total_peaks = 0
    total_peaks_pos = []

    ipeaks_int = peak_local_max(colcut, min_distance=2, threshold_rel=0.2)[:, 0]
    ipeaks_float = refine_peaks(colcut, ipeaks_int, 3)[0]
    peaks_y = np.ones((ipeaks_int.shape[0],3))
    peaks_y[:,0] = ipeaks_int
    peaks_y[:,1] = ipeaks_float
    peaks_y[:,2] = colcut[ipeaks_int]
    box_match = np.digitize(peaks_y[:, 0], box_borders)

    _logger.debug('pairing fibers')
    for box in boxes:
        nfibers = box['nfibers']
        boxid = box['id'] - 1

        dist_b_fibs = (box_borders[boxid + 1] - box_borders[boxid]) / (nfibers + 2.0)
        mask_fibers = (box_match == (boxid + 1))
        # Peaks in this box
        thispeaks = peaks_y[mask_fibers]
        npeaks = len(thispeaks)
        total_peaks += npeaks
        for elem in thispeaks:
            total_peaks_pos.append(elem.tolist())

        _logger.debug('box: %s', box['id'])
        # Start by matching the first peak
        # with the first fiber
        fid = 0
        current_peak = 0
        pairs_1 = [(fid, current_peak)]
        fid += 1

        scale = 1
        while (current_peak < npeaks - 1) and (fid < nfibers):
            # Expected distance to next fiber
            expected_distance = scale * dist_b_fibs
            _logger.debug('expected %s', expected_distance)
            _logger.debug('current peak: %s', current_peak)
            for idx in range(current_peak + 1, npeaks):
                distance = abs(thispeaks[idx, 1] - thispeaks[current_peak, 1])
                if abs(distance - expected_distance) <= tol:
                    # We have a match
                    # We could update
                    # dist_b_fibs = distance / scale
                    # But is not clear this is better

                    # Store this match
                    pairs_1.append((fid, idx))
                    current_peak = idx
                    # Next
                    scale = 1
                    break
            else:
                # This fiber has no match
                pairs_1.append((fid, None))
                # Try a fiber further away
                scale += 1
            # Match next fiber
            fid += 1
        _logger.debug(pairs_1)
        _logger.debug('matched %s \t missing: %s', len(pairs_1),nfibers-len(pairs_1))
        remainig = nfibers - len(pairs_1)
        if remainig > 0:
            _logger.debug('We have to pair: %s', remainig)
            # Position of first match fiber

            # Position of last match fiber
            for fid, peakid in reversed(pairs_1):
                if peakid is not None:
                    last_matched_peak = peakid
                    last_matched_fiber = fid
                    break
            else:
                raise ValueError('None matched')
            _logger.debug('peaks: %s \t %s', thispeaks[0, 1], thispeaks[last_matched_peak, 1])
            _logger.debug('borders: %s \t %s', box_borders[boxid], box_borders[boxid+1])
            ldist = thispeaks[0, 1] - box_borders[boxid]
            rdist = box_borders[boxid + 1] - thispeaks[last_matched_peak, 1]
            lcap = ldist / dist_b_fibs - 1
            rcap = rdist / dist_b_fibs - 1
            _logger.debug('L distance %s \t %s', ldist, lcap)
            _logger.debug('R distance %s \t %s', rdist, rcap)
            lcapi = int(lcap + 0.5)
            rcapi = int(rcap + 0.5)

            on_r = rcapi <= lcapi
            mincap = min(lcapi, rcapi)
            maxcap = max(lcapi, rcapi)

            cap1 = min(mincap, remainig)
            cap2 = min(maxcap, remainig - cap1)
            cap3 = remainig - cap1 - cap2

            if cap3 > 0:
                _logger.debug('we dont have space %s fibers no allocated', cap3)

            if on_r:
                # Fill rcap fibers, then lcap
                capr = cap1
                capl = cap2
            else:
                capr = cap2
                capl = cap1

            addl = [(x, None) for x in range(-capl, 0)]
            addr = [(x, None) for x in range(last_matched_fiber + 1, last_matched_fiber + 1 + capr)]
            _logger.debug('add %s fibers on the right', capr)
            _logger.debug('add %s fibers on the left', capl)
            _logger.debug(addl)
            _logger.debug(addr)
            pairs_1 = addl + pairs_1 + addr
            _logger.debug(pairs_1)

        # reindex
        # assert(len(pairs_1) == nfibers)

        for fibid, (relfibid, match) in enumerate(pairs_1, counted_fibers):
            fiber_traces[fibid] = FiberTraceInfo(fibid+1, box['id'])
            if match is not None:
                fiber_traces[fibid].start = (center, thispeaks[match, 1], thispeaks[match, 2])
            # else:
            #     fiber_traces[fibid].start = (center, 0, 0)
        counted_fibers += nfibers

        # import matplotlib.pyplot as plt
        # plt.xlim([box_borders[boxid], box_borders[boxid + 1]])
        # plt.plot(colcut, 'b-')
        # plt.plot(thispeaks[:, 1], thispeaks[:, 2], 'ro')
        # plt.plot(peaks_y[:,1], peaks_y[:, 2], 'ro')
        # plt.title('Box %s' %box['id'])
        # plt.show()

    # import matplotlib.pyplot as plt
    # total_peaks_pos = np.array(total_peaks_pos)
    # plt.plot(colcut, 'b-')
    # plt.plot(total_peaks_pos[:, 1], total_peaks_pos[:, 2], 'ro')
    # plt.show()

    _logger.debug ('total found peaks: %s' %total_peaks)
    _logger.debug ('total found + recovered peaks: %s' %counted_fibers)

    return fiber_traces
