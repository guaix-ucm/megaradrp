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

from .peakdetection import peak_detection_mean_window  #Should be avoided
from numina.array.peaks.peakdet import find_peaks_indexes, refine_peaks

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

    box_borders = [123, 262, 400, 536, 713, 887, 1103, 1360, 1823, 2287, 2750, 3007, 3223, 3397, 3574, 3710, 3848, 3985]
    box_borders2 = [150, 294, 429, 557, 733, 906, 1120, 1373, 1823, 2287, 2737, 2992, 3205, 3377, 3553, 3684, 3815, 3962]
    vborders2 = [676, 818, 1966, 1222, 3502, 5415, 2110, 2850, 2133, 2800, 3004, 4134, 6046, 4886, 1458, 2614, 900, 721]

    for x1, x2 in zip(box_borders, box_borders2):
        print (x1, x1-x2)

    xx = np.arange(image.shape[0])

    cut_region = slice(center-hs, center+hs)
    cut = image[:,cut_region]
    colcut = cut.mean(axis=1)

    import matplotlib.pyplot as plt

    print('background', background)
    maxt = peak_detection_mean_window(colcut, x=xx, k=3, xmin=ixmin, xmax=ixmax, background=background)
    print('maxt', maxt)
    #plt.ylim([-600, 600])
    plt.plot(colcut)
    #plt.plot(maxt[:,0], maxt[:,2], 'r*')
    plt.plot(box_borders, [5000 for i in box_borders], 'g*')
    plt.plot(box_borders2, vborders2, 'r*')
    plt.show()
    peakdist = np.diff(maxt[:,1])
    npeaks = maxt.shape[0]

    # Count the number of groups
    # fibers separated by more than "maxdis" pixels
    group_changes = np.nonzero(peakdist > maxdis)[0]
    # Fill groups
    groups = np.ones(npeaks, dtype='int')
    for p in group_changes:
        groups[p+1:] += 1

    fiber_traces = {}
    for fibid in range(1, npeaks+1):
        fiber_traces[fibid] = FiberTraceInfo(fibid, int(groups[fibid-1]))

    # Window to fix the center of the trace
    fw = 2
    for fibid, trace in fiber_traces.items():
        _xi, _yi, _vali = center, maxt[fibid-1,1], maxt[fibid-1,2]
        pixmax = int(maxt[fibid-1,0])
        # Take 2*2+1 pix
        # This part and interp_max_3(image[nearp3-1:nearp3+2, col])
        # should do the same
        tx, py, _pos = delicate_centre(xx[pixmax-fw: pixmax+fw+1], 
                                     colcut[pixmax-fw: pixmax+fw+1])

        trace.start = (center, tx, py)

    return fiber_traces


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

def ncl_gen_triplets_master(positions):
    import itertools

    nlines_master = len(positions)
    iter_comb_triplets = itertools.combinations(range(nlines_master), 3)
    triplets_master_list = list(iter_comb_triplets)
    # print(iter_comb_triplets) This is iterator
    # For each triplet, compute the relative position of the central line.

    ntriplets_master = len(triplets_master_list)

    ratios_master = np.zeros(ntriplets_master)
    for i_tupla in range(ntriplets_master):
        i1, i2, i3 = triplets_master_list[i_tupla]
        delta1 = positions[i2] - positions[i1]
        delta2 = positions[i3] - positions[i1]
        ratios_master[i_tupla] = delta1 / delta2

    # Compute the array of indices that index the above ratios in sorted order.
    isort_ratios_master = np.argsort(ratios_master)

    # Simultaneous sort of position ratios and triplets.
    ratios_master_sorted = ratios_master[isort_ratios_master]
    triplets_master_sorted_list = [triplets_master_list[i] for i in isort_ratios_master]

    return ntriplets_master, ratios_master_sorted, triplets_master_sorted_list




def estimate_thresholds(image, center, hs, boxref):
    """Estimate background from values in boxes between fibers"""

    cut_region = slice(center-hs, center+hs)
    cut = image[boxref, cut_region]

    colcut = cut.mean(axis=1)
    colcut_std = cut.std(axis=1)

    max_val = colcut
    max_std = colcut_std
    background = max_val + 2 * max_std

    return background

def init_traces_ex(image, center, hs, box_borders, tol=1.5):

    cut_region = slice(center-hs, center+hs)
    cut = image[:,cut_region]
    colcut = cut.mean(axis=1)

    ##########################################################
    # Iba aqui
    ##########################################################

    counted_fibers = 0
    fiber_traces = {}
    total_peaks = 0
    total_peaks_pos = []

    thresholds = estimate_thresholds(image, center, 3, box_borders)

    # print('pairing fibers')
    for box in boxes:
        nfibers = box['nfibers']
        boxid = box['id'] - 1

        ##########################################################

        ipeaks_int = find_peaks_indexes(colcut, 3, thresholds[boxid])
        ipeaks_float = refine_peaks(colcut, ipeaks_int, 3)[0]
        peaks_y = np.ones((ipeaks_int.shape[0],3))
        peaks_y[:,0] = ipeaks_int
        peaks_y[:,1] = ipeaks_float
        peaks_y[:,2] = colcut[ipeaks_int]


        # Interval where the peak is
        box_match = np.digitize(peaks_y[:, 0], box_borders)
        ##########################################################


        dist_b_fibs = (box_borders[boxid + 1] - box_borders[boxid]) / (nfibers + 2.0)
        mask_fibers = (box_match == (boxid + 1))
        # Peaks in this box
        thispeaks = peaks_y[mask_fibers]
        npeaks = len(thispeaks)
        total_peaks += npeaks
        for elem in thispeaks:
            total_peaks_pos.append(elem.tolist())



        plt_expected_pos = box_borders[boxid] + dist_b_fibs* np.arange(0, nfibers + 2)

        print('box:', box['id'])
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
            print('expected ', expected_distance)
            print('current peak', current_peak)
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
        print(pairs_1)
        print('matched', len(pairs_1), 'missing', nfibers-len(pairs_1))
        remainig = nfibers - len(pairs_1)
        if remainig > 0:
            print('We have to pair', remainig)
            # Position of first match fiber

            # Position of last match fiber
            for fid, peakid in reversed(pairs_1):
                if peakid is not None:
                    last_matched_peak = peakid
                    last_matched_fiber = fid
                    break
            else:
                raise ValueError('None matched')
            print('peaks', thispeaks[0, 1], thispeaks[last_matched_peak, 1])
            print('borders', box_borders[boxid], box_borders[boxid+1])
            ldist = thispeaks[0, 1] - box_borders[boxid]
            rdist = box_borders[boxid + 1] - thispeaks[last_matched_peak, 1]
            lcap = ldist / dist_b_fibs - 1
            rcap = rdist / dist_b_fibs - 1
            print('L distance', ldist, lcap)
            print('R distance', rdist, rcap)
            lcapi = int(lcap + 0.5)
            rcapi = int(rcap + 0.5)

            on_r = rcapi <= lcapi
            mincap = min(lcapi, rcapi)
            maxcap = max(lcapi, rcapi)

            cap1 = min(mincap, remainig)
            cap2 = min(maxcap, remainig - cap1)
            cap3 = remainig - cap1 - cap2

            if cap3 > 0:
                print('we dont have space', cap3, 'fibers no allocated')

            if on_r:
                # Fill rcap fibers, then lcap
                capr = cap1
                capl = cap2
            else:
                capr = cap2
                capl = cap1

            addl = [(x, None) for x in range(-capl, 0)]
            addr = [(x, None) for x in range(last_matched_fiber + 1, last_matched_fiber + 1 + capr)]
            print('add', capr, 'fibers on the right')
            print('add', capl, 'fibers on the left')
            print(addl)
            print(addr)
            pairs_1 = addl + pairs_1 + addr
            print(pairs_1)

        # reindex
        assert(len(pairs_1) == nfibers)

        for fibid, (relfibid, match) in enumerate(pairs_1, counted_fibers):
            if match is not None:
                fiber_traces[fibid] = FiberTraceInfo(fibid, box['id'])
                fiber_traces[fibid].start = (center, thispeaks[match, 1], thispeaks[match, 2])
            else:
                # FIXME: what to do for no matches
                pass
        counted_fibers += nfibers


        # import matplotlib.pyplot as plt
        # plt.xlim([box_borders[boxid], box_borders[boxid + 1]])
        # plt.plot(colcut, 'b-')
        # plt.vlines(plt_expected_pos[1:-1],0,40000)
        # plt.plot(thispeaks[:, 1], thispeaks[:, 2], 'ro')
        # plt.plot(peaks_y[:,1], peaks_y[:, 2], 'ro')
        # plt.plot(plt_expected_pos, [(1.1 * thispeaks[:, 2].max()) for _ in plt_expected_pos], 'go')
        # plt.title('Box %s' %box['id'])
        # plt.show()

        print ('*' * 250)

    print ('total found peaks: %s' %total_peaks)

    import matplotlib.pyplot as plt
    total_peaks_pos = np.array(total_peaks_pos)
    plt.plot(colcut, 'b-')
    plt.plot(total_peaks_pos[:, 1], total_peaks_pos[:, 2], 'ro')
    plt.show()

    return fiber_traces
