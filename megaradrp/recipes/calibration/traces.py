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

from megaradrp.trace.peakdetection import peak_detection_mean_window

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
    maxt = peak_detection_mean_window(y, x=x, k=3,xmin=ixmin,xmax=ixmax,background=background)
    npeaks = len(maxt)
    peakdist = np.diff(maxt[:, 1])
    # number of peaks
    _logger.debug('npeaks %d', npeaks)
    _logger.debug('distances between peaks (max) % f (min) %f', peakdist.max(), peakdist.min())
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
        tx, py, _pos = delicate_centre(x[pixmax-2: pixmax+2+1], y[pixmax-2: pixmax+2+1])
        trace.sample_f.append(center)
        trace.trace_f.append(tx)
        trace.peak_f.append(py)
    return fiber_traces


