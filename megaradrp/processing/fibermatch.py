#
# Copyright 2011-2018 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#


"""Match and identify fibers"""

import itertools

from six.moves import filter as ifilter


FIBER_PEAK, FIBER_DEAD = (0, 1)
FOUND_PEAK, FOUND_VALLEY, EXPECT_VALLEY = (0, 1, 2)


def generate_box_model(nfibers, startid=1, scale=6.0, start=0.0,
                   missing=None, skipid=None,
                   missing_count='rel'):
    """Generate a model of the expected peaks in a box"""

    if missing_count == 'rel':
        off = 1
    elif missing_count == 'abs':
        off = 0
    else:
        raise ValueError('missing_count must be either "rel" or "abs"')

    if skipid is None:
        skipid = []

    if missing is None:
        missing = []

    iter1 = itertools.count(startid)
    iter2 = ifilter(lambda x:x not in skipid, iter1)
    iter3 = itertools.islice(iter2, nfibers)

    result = []
    for idx, val in zip(iter3, range(nfibers)):
        key = FIBER_PEAK
        if off and (val + 1) in missing:
            key = FIBER_DEAD
        result.append((idx, val, start + scale * val , key))
    return result


def count_peaks(peaks, tol=1.2, distance=6.0, start=1):
    """Count the peaks and valleys of an array"""

    expected_distance = distance
    scale = 1
    max_scale_jum = 3

    pidx = 0
    pid = start

    if len(peaks) == 0:
        raise ValueError('no peaks to count')

    p1, rest = peaks[0], peaks[1:]
    # pref = p1
    values = [(pid, p1, FOUND_PEAK, 0)]
    measured_dists = []

    while len(rest) > 0:
        # print('im peak:', p1, ' next peak should be around:', p1 + scale * expected_distance)
        p2 = rest[0]
        dist = abs(p1 - p2)
        while True:
            sed = scale * expected_distance
            # print('next peak is:', p2, 'distance from p1 is', dist)
            # print('expected distance is:', sed)
            pid += 1
            pidx += 1
            if abs(dist - sed) < scale * tol:
                # print('p2 is within expected distance with tol:', tol)
                # print('p2 is next peak, with scale', scale)
                new_exp = dist / scale
                if abs(new_exp - expected_distance) < 1:
                    expected_distance = dist / scale

                measured_dists.append(expected_distance)
                values.append((pid, p2, FOUND_PEAK, pidx))
                scale = 1
                break
            else:
                # print('p2 is not within expected distance with tol:', tol)
                pex = p1 + sed
                values.append((pid, pex, FOUND_VALLEY, None))
                scale += 1
                # print('increase scale to:', scale)
                if scale > max_scale_jum:
                    # print('moving to far')
                    raise ValueError('moving too far apart')

        p1, rest = rest[0], rest[1:]
    return values, measured_dists


# matching
def valid_match(model, values):
    """Match a model with found peaks"""
    for m, v in zip(model, values):
        field_m = m[3]
        field_v = v[2]
        if field_m == FIBER_DEAD and field_v == FOUND_PEAK:
            # print(m, v, 'must be empty')
            return False
    return True


def complete_solutions(model, values, borders, scale=6.0):
    """Try different solutions to match the missing peaks"""

    border1 = borders[0]
    border2 = borders[1]

    missing = len(model) - len(values)

    # last_peak = values[-1][1]
    last_peak = values[-1]
    first_peak = values[0]
    # distance to borders
    d1 = abs(first_peak[1] - border1)
    d2 = abs(last_peak[1] - border2)

    solutions = []

    # print('XXX missing peaks', missing)
    for pl in range(missing + 1):
        pr = missing - pl
        # print('XXXX', pl, pr)
        # Complete peaks
        # preapend pl peaks

        pre = []
        post = []
        # print('XXXX space per peak added', d1 / (pl + 1.0), d2 / (pr + 1.0))
        for peak in range(pl):
            idx = -peak - 1
            dist = scale * idx
            tok = (idx, dist, EXPECT_VALLEY, None)
            pre.append(tok)

        for peak in range(pr):
            idx = (len(values) + peak) + 1
            dist = last_peak[1] + scale * (peak + 1)
            tok = (idx, dist, EXPECT_VALLEY, None)
            post.append(tok)

        if valid_match(model, itertools.chain(pre, values, post)):
            # print('computing score')
            # print(border1, border2)
            ratios = [d1 / (pl +1 ), d2 / (pr + 1)]
            # print('ratios:', ratios)
            exp_dist = scale
            score = sum([(c - exp_dist) ** 2 for c in ratios])
            # print('score', score)
            solution = (score, (pre, post))
            solutions.append(solution)

    # print('------------')
    return solutions


def iter_best_solution(model, values, solutions):

    # This works because score is the first...
    best_solution = min(solutions)
    minscore, added_peaks = best_solution

    pre, post = added_peaks
    comp_sol = itertools.chain(pre, values, post)
    for xx, yy in zip(comp_sol, model):
        fibid = yy[0]
        match = xx[3]
        yield (fibid, match)