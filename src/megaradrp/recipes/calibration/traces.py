#
# Copyright 2014-2015 Universidad Complutense de Madrid
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

'''Peak tracing of spectral profiles'''

from __future__ import division

import logging

import numpy as np
from scipy.interpolate import UnivariateSpline

from .peakdetection import peak_detection_mean_window

_logger = logging.getLogger('megara.trace')


def delicate_centre(x, y):
    pos = np.polyfit(x, y, deg=2)
    tx = -pos[1] / (2 * pos[0])
    py = np.polyval(pos, tx)
    return tx, py, pos


class Trace(object):
    def __init__(self, start=0):
        self.start = start
        self.lost = None
        self.sample_c = []
        self.trace_c = []
        self.peak_c = []
        self.sample_f = []
        self.trace_f = []
        self.peak_f = []

    def predict_position(self, col):
        return self.trace_f[-1]


class FiberTrace(Trace):
    def __init__(self, fibid, start=0):
        super(FiberTrace, self).__init__(start=start)
        self.fibid = fibid


def init_traces(y, center, x, ixmin, ixmax, background):
    maxt = peak_detection_mean_window(y, x=x, k=3,
                                      xmin=ixmin,
                                      xmax=ixmax,
                                      background=background)
    npeaks = len(maxt)
    peakdist = np.diff(maxt[:, 1])
    # number of peaks
    _logger.debug('npeaks %d', npeaks)
    _logger.debug('distances between peaks (max) % f (min) %f',
                  peakdist.max(), peakdist.min())
    #
    fiber_traces = {i: FiberTrace(i, start=center) for i in range(npeaks)}
    #
    # peaks

    for fibid, trace in fiber_traces.items():
        trace.sample_c.append(center)
        trace.trace_c.append(maxt[fibid, 1])
        trace.peak_c.append(maxt[fibid, 2])
        pixmax = int(maxt[fibid, 0])
        # Take 2*2+1 pix
        tx, py, _pos = delicate_centre(x[pixmax-2: pixmax+2+1],
                                       y[pixmax-2: pixmax+2+1])
        trace.sample_f.append(center)
        trace.trace_f.append(tx)
        trace.peak_f.append(py)
    return fiber_traces


def explore_traces_left(data, x, cstart, cend, hs, fiber_traces,
                        background, maxdis):
    return explore_traces_common(data, x, cstart, cend, -hs, fiber_traces,
                                 background, maxdis)


def explore_traces_rigth(data, x, cstart, cend, hs, fiber_traces,
                         background, maxdis):
    return explore_traces_common(data, x, cstart, cend, hs, fiber_traces,
                                 background, maxdis)


def explore_traces_common(data, x, cstart, cend, hs, fiber_traces,
                          background, maxdis):

    centers = range(cstart, cend, hs)

    hs = abs(hs)
    for c in centers[1:]:
        lost_count = 0

        cut_region = slice(c-hs, c+hs)
        colcut = data[:, cut_region].mean(axis=1)
        # FIXME: background could be variable...
        maxt = peak_detection_mean_window(colcut, x=x, k=3,
                                          background=background)
        npeaks = len(maxt)
        _logger.debug('npeaks =%i, c=%i, hs=%i', npeaks, c, hs)
        if npeaks == 0:
            _logger.debug('no more peaks')
            for fibid, trace in fiber_traces.items():
                if trace.lost is None:
                    trace.lost = c
            break

        # pixel fraction peak
        # fit peaks
        fitt = maxt.copy()
        for peakid in range(npeaks):
            pixmax = int(maxt[peakid, 0])
            # Take 2*2+1 pix
            xi = x[pixmax-2: pixmax+2+1]
            yi = colcut[pixmax-2: pixmax+2+1]
            if xi.size < 5 or yi.size < 5:
                fitt[peakid] = maxt[peakid]
            else:
                tx, py, pos = delicate_centre(xi, yi)
                fitt[peakid, 1:] = [tx, py]
        start = 0
        for fibid, trace in fiber_traces.items():
            if trace.lost is not None:
                lost_count += 1
                continue
            _logger.debug('prev fib=%i is peak=%i', fibid, peakid)

            # predict position for this trace in this position
            pos = trace.predict_position(c)

            _logger.debug('predicted position of fib=%i is pos=%f', fibid, pos)
            _logger.debug('find current peak %f pix near pos=%f', maxdis, pos)
            prev_dis = np.inf

            for peakid in range(start, npeaks):
                peak_pos = fitt[peakid, 1]
                peak_val = fitt[peakid, 2]
                dis = abs(peak_pos - pos)
                if dis > prev_dis:
                    # We are not getting closer, break here
                    # trace is lost
                    _logger.debug('peak %i is farther than initial '
                                  'peak %i from fib=%i', peakid, start, fibid)
                    lost_count += 1
                    trace.lost = c
                    break
                else:
                    prev_dis = dis

                if dis < maxdis:
                    _logger.debug('peak %i pos %f is within %f '
                                  'pixels of pos=%f fib=%i',
                                  peakid, peak_pos, maxdis, peak_pos, fibid)
                    trace.sample_c.append(c)
                    trace.sample_f.append(c)
                    trace.trace_f.append(peak_pos)
                    trace.peak_f.append(peak_val)
                    trace.trace_c.append(maxt[peakid, 1])
                    trace.peak_c.append(maxt[peakid, 2])
                    start = peakid + 1
                    break
                else:
                    pass
                    _logger.debug('peak %i is not within %f '
                                  'pixels of pos=%f fib=%i',
                                  peakid, maxdis, peak_pos, fibid)
            else:
                _logger.debug('all peaks are farther than '
                              'initial peak %i from fib=%i', start, fibid)
                lost_count += 1
                trace.lost = c

        _logger.debug('col = %i peaks = %i lost fibers = %i',
                      c, npeaks, lost_count)


def domefun(image, maxdis=1.5, hs=5, cstart=None):
    '''
        hs: half size of the cut region
    '''

    wl_axis = 1
    sp_axis = 0

    # Number of steps needed to predict the position
    # npredict = 0

    shape = image.shape
    wl_len = shape[wl_axis]
    sp_len = shape[sp_axis]

    xmin = 0
    xmax = sp_len

    if cstart is None:
        cstart = (xmax - xmin) // 2

    data = image[xmin:xmax]

    # Use this to avoid bad pixels in the borders
    # No peak will be detected outside this
    ixmin = 10
    ixmax = sp_len - 10

    # No peak will be detected below this flux
    # TODO: change this
    background = 10.0

    # Init traces
    sp_x = np.arange(xmin, xmax)

    cut_region = slice(cstart-hs, cstart+hs)
    colcut = data[:, cut_region].mean(axis=1)

    fiber_traces = init_traces(colcut, cstart, sp_x,
                               ixmin, ixmax, background)

    cend1 = 50
    cend2 = wl_len

    explore_traces_left(data, sp_x, cstart, cend1, hs,
                        fiber_traces, background, maxdis)

    def reverse_trace(trace):
        trace.lost = None
        trace.sample_c.reverse()
        trace.sample_f.reverse()
        trace.trace_c.reverse()
        trace.trace_f.reverse()
        trace.peak_c.reverse()
        trace.peak_f.reverse()
        return trace

    for trace in fiber_traces.values():
        reverse_trace(trace)

    explore_traces_rigth(data, sp_x, cstart, cend2, hs,
                         fiber_traces, background, maxdis)

    fit_traces = {key: fit_trace(trace, 5)
                  for (key, trace) in fiber_traces.items()}

    return fit_traces


def fit_trace(trace, deg=5):
    spl = UnivariateSpline(trace.sample_f, trace.trace_f, k=deg)
    return [deg, spl.get_coeffs().tolist(), spl.get_knots().tolist()]
