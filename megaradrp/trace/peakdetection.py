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

'''Peak finding for Megara'''

import numpy as np
from scipy.ndimage.filters import generic_filter

def _vecS2(k, data):
    '''min filter for peak detection.'''

    def func(x):
        return x[k] - 0.5 * (x[:k].mean() + x[k+1:].mean())

    ap = generic_filter(data, func, size=2*k+1)
    return ap


# http://tcs-trddc.com/trddc_website/pdf/srl/palshikar_sapdts_2009.pdf
def _generic_peak_filtering(method, y, x=None, k=3, h=1.0, background=0.0,xmin=None, xmax=None):
    '''Helper class for filtered peak finding.'''

    if x is None:
        x = np.arange(len(y))

    if xmin is None:
        xmin = x[0]

    if xmax is None:
        xmax = x[-1]

    candidates = []
    # Fixme...

    ap = method(k, y)
    # mean and std of pos values
    # ppos = ap[ap>0]
    #
    # This does not work, we need something local, not global
    # mpos = appos.mean()
    # spos = appos.std()

    for idx, a in enumerate(ap):
        if a > 0 and y[idx] > background and (xmin <= x[idx] <= xmax):
            candidates.append(idx)

    # Filter close elements
    result = []
    while len(candidates) > 1:

        x0 = candidates[0]
        for idx, xi in enumerate(candidates[1:], 1):
            if abs(xi - x0) > k:
                break

        group = candidates[:idx]

        remaining = candidates[idx:]
        maxvalid = y[group].argmax()
        finalid = group[maxvalid]
        result.append((finalid, x[finalid], y[finalid]))

        candidates = remaining

    return np.array(result)


def peak_detection_mean_window(y, x=None, k=3, h=1.0, background=0.0, xmin=None, xmax=None):
    '''Detect peaks using a mean filter in a 2*k +1 window.

    The filter value is the central point minus the mean of the 2*k
    adjacent points.
    '''
    return _generic_peak_filtering(_vecS2, y, x=x, k=k, h=h,background=background,xmin=xmin, xmax=xmax)



def peakdet(v, delta, x=None, back=0.0):
    '''Basic peak detection.'''

    if x is None:
        x = np.arange(len(v))

    v = np.asarray(v)

    if len(v) != len(x):
        raise ValueError('Input vectors v and x must have same length')

    if not np.isscalar(delta):
        raise TypeError('Input argument delta must be a scalar')

    if delta <= 0:
        raise ValueError('Input argument delta must be positive')

    maxtab = []
    mintab = []

    mn, mx = np.Inf, -np.Inf
    mnpos, mxpos = np.NaN, np.NaN

    lookformax = True

    for i, this in enumerate(v):
        if this > mx and this > back:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]

        if lookformax:
            if this < mx - delta and this > back:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn + delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True

    return np.array(maxtab), np.array(mintab)