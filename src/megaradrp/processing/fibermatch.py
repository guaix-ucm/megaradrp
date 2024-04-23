#
# Copyright 2011-2024 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#


"""Match and identify fibers"""
import enum
import itertools

import attrs


class PeakMode(enum.Enum):
    FIBER_PEAK = 0
    FIBER_DEAD = 1


class PeakFound(enum.Enum):
    FOUND_PEAK = 0
    FOUND_VALLEY = 1
    EXPECT_VALLEY = 2
    MAX_JUMP = 3


@attrs.define
class FiberModelElement:
    fibid: int = attrs.field()
    mode: PeakMode = attrs.field()


@attrs.define
class PeakMatch:
    count: int = attrs.field()
    pos: float = attrs.field()
    mode: PeakFound = attrs.field()
    idx: int = attrs.field()


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
    iter2 = filter(lambda x: x not in skip_fibids, iter1)
    iter3 = itertools.islice(iter2, nfibers)

    result = []
    for idx, fibid in enumerate(iter3, 1):
        key = PeakMode.FIBER_PEAK
        if idx in missing_relids:
            key = PeakMode.FIBER_DEAD
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

    values = [PeakMatch(pid, p1, PeakFound.FOUND_PEAK, 0)]

    while len(rest) > 0:
        # print(f'im peak: {p1:.2f}, next peak should be around: {p1 + scale * expected_distance:.2f}')
        p2 = rest[0]
        dist = abs(p1 - p2)
        last_info = values[-1]
        while True:
            sed = scale * expected_distance
            # print(f'next peak is: {p2:.2f}, distance from p1 is {dist:.2f}')
            # print('expected distance is:', sed)
            pid += 1

            if abs(dist - sed) < scale * tol:
                # print(f'p2 is within expected distance with tol {tol}')
                # print(f'p2 is next peak, with scale', scale)
                pidx += 1
                values.append(PeakMatch(pid, p2, PeakFound.FOUND_PEAK, pidx))
                scale = 1

                p1, rest = rest[0], rest[1:]
                break
            else:
                # print('p2 is not within expected distance with tol:', tol)
                pex = p1 + sed
                pidx += 1
                values.append(PeakMatch(pid, pex, PeakFound.FOUND_VALLEY, None))
                scale += 1
                # print('increase scale to:', scale)
                # Skip this peak
                rest = rest[1:]
                if scale > max_scale_jump:
                    # print('moving too far away')
                    msg = f'peak {pid} not found within expected distance from peak {last_info.count}'
                    print(msg)
                    # end process
                    rest = rest[-1:]
                break
    return values


# matching
def valid_match(model, values):
    """Match a model with found peaks"""
    for m, v in zip(model, values):
        field_m = m.mode
        field_v = v.mode
        if field_m == PeakMode.FIBER_DEAD and field_v == PeakFound.FOUND_PEAK:
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
    d1 = abs(first_peak.pos - border1)
    d2 = abs(last_peak.pos - border2)

    solutions = []

    # print('XXX missing peaks', missing)
    for pl in range(missing + 1):
        pr = missing - pl
        # print('XXXX', pl, pr)
        # Complete peaks
        # preappend pl peaks

        pre = []
        post = []
        # print('XXXX space per peak added', d1 / (pl + 1.0), d2 / (pr + 1.0))
        for peak in range(pl):
            idx = -peak - 1
            dist = scale * idx
            tok = PeakMatch(idx, dist, PeakFound.EXPECT_VALLEY, None)
            pre.append(tok)

        for peak in range(pr):
            idx = (len(values) + peak) + 1
            dist = last_peak.pos + scale * (peak + 1)
            tok = PeakMatch(idx, dist, PeakFound.EXPECT_VALLEY, None)
            post.append(tok)

        if valid_match(model, itertools.chain(pre, values, post)):
            # print('computing score')
            # print(border1, border2)
            ratios = [d1 / (pl + 1), d2 / (pr + 1)]
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
        fibid = yy.fibid
        match = xx.idx
        yield fibid, match
