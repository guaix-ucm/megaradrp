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

from .peakdetection import peak_detection_mean_window

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


def init_traces(image, center, hs, background, npred, maxdis=9.0):

    ixmin= 0
    ixmax = image.shape[0]
    
    xx = np.arange(image.shape[0])

    cut_region = slice(center-hs, center+hs)
    cut = image[:,cut_region]
    colcut = cut.mean(axis=1)
    maxt = peak_detection_mean_window(colcut, x=xx, k=3, xmin=ixmin, xmax=ixmax, background=background)
    #npeaks = len(maxt)
    peakdist = np.diff(maxt[:,1])
    # number of peaks    #
    fiber_traces = {}
    fibid = 1
    gcounter = 0
    boxid = 1
    for _xpeak, dis in zip(maxt[:,1], peakdist):
        gcounter += 1
        if dis > maxdis:
            #print(xpeak, dis, gcounter, boxid, fibid)
            gcounter = 0
            boxid += 1
        fiber_traces[fibid] = FiberTraceInfo(fibid, boxid)
        fibid += 1

    fw = 2

    for fibid, trace in fiber_traces.items():
        _xi, _yi, _vali = center, maxt[fibid,1], maxt[fibid,2]
        pixmax = int(maxt[fibid,0])
        # Take 2*2+1 pix
        # This part and interp_max_3(image[nearp3-1:nearp3+2, col])
        # should do the same
        tx, py, _pos = delicate_centre(xx[pixmax-fw: pixmax+fw+1], 
                                     colcut[pixmax-fw: pixmax+fw+1])

        trace.start = (center, tx, py)

    return fiber_traces
