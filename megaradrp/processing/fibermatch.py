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
import collections


from six.moves import filter as ifilter


FIBER_PEAK, FIBER_DEAD = (0, 1)
FOUND_PEAK, FOUND_VALLEY, EXPECT_VALLEY = (0, 1, 2)


FiberModelElement = collections.namedtuple('FiberModelElement', ['fibid', 'mode'])


def generate_box_model(nfibers, start=1,
                       missing_relids=None,
                       skip_fibids=None
                       ):
    """Generate a model of the expected peaks in a box"""

    if skip_fibids is None:
        skip_fibids = []

    if missing_relids is None:
        missing_relids = []

    iter1 = itertools.count(start)
    iter2 = ifilter(lambda x:x not in skip_fibids, iter1)
    iter3 = itertools.islice(iter2, nfibers)

    result = []
    for idx, fibid in enumerate(iter3, 1):
        key = FIBER_PEAK
        if idx in missing_relids:
            key = FIBER_DEAD
        tok = FiberModelElement(fibid=fibid, mode=key)
        result.append(tok)
    return result


def count_peaks(peaks, tol=1.2, distance=6.0, start=1, max_scale_jump=3):
    """Count the peaks and valleys of an array"""

    expected_distance = distance

    scale = 1
    pidx = 0
    pid = start

    if len(peaks) == 0:
        raise ValueError('no peaks to count')

    p1, rest = peaks[0], peaks[1:]
    # pref = p1
    values = [(pid, p1, FOUND_PEAK, 0)]

    while len(rest) > 0:
        # print('im peak:', p1, ' next peak should be around:', p1 + scale * expected_distance)
        p2 = rest[0]
        dist = abs(p1 - p2)
        while True:
            sed = scale * expected_distance
            # print('next peak is:', p2, 'distance from p1 is', dist)
            # print('expected distance is:', sed)
            pid += 1

            if abs(dist - sed) < scale * tol:
                # print('p2 is within expected distance with tol:', tol)
                # print('p2 is next peak, with scale', scale)
                pidx += 1
                values.append((pid, p2, FOUND_PEAK, pidx))
                scale = 1
                break
            else:
                # print('p2 is not within expected distance with tol:', tol)
                pex = p1 + sed
                values.append((pid, pex, FOUND_VALLEY, None))
                scale += 1
                # print('increase scale to:', scale)
                if scale > max_scale_jump:
                    # print('moving to far')
                    raise ValueError('moving too far apart')

        p1, rest = rest[0], rest[1:]
    return values


# matching
def valid_match(model, values):
    """Match a model with found peaks"""
    for m, v in zip(model, values):
        field_m = m.mode
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