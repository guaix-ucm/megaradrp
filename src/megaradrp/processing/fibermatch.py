#
# Copyright 2011-2024 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#


"""Match and identify fibers"""

import logging
import enum
import itertools

import attrs

_logger = logging.getLogger(__name__)


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


def generate_box_model(nfibers, start=1, missing_relids=None, skip_fibids=None):
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
    """Count the peaks and valleys of an array

    The peaks are expected every `distance` pixels. A peak at k times
    `distance` from the previous one (within k times `tol`) leaves k - 1
    valleys, missing peaks, before it. A peak closer than expected to the
    previous one is skipped. The counting stops when the next peak is
    `max_scale_jump` or more times `distance` away.

    Parameters
    ----------
    peaks : array_like
        Sorted positions of the peaks
    tol : float
        Tolerance in the distance between consecutive peaks
    distance : float
        Expected distance between consecutive peaks
    start : int
        Count of the first peak
    max_scale_jump : int
        Distance, in units of `distance`, at which the counting stops

    Returns
    -------
    list of PeakMatch
        The peaks and the valleys between them. The `idx` of a peak
        is its index in `peaks`.
    """

    if len(peaks) == 0:
        raise ValueError("no peaks to count")

    pid = start
    values = [PeakMatch(pid, peaks[0], PeakFound.FOUND_PEAK, 0)]
    valleys = []
    p1 = peaks[0]
    scale = 1
    idx = 1
    while idx < len(peaks):
        p2 = peaks[idx]
        dist = p2 - p1
        sed = scale * distance
        if abs(dist - sed) < scale * tol:
            # the next peak, after scale - 1 valleys
            values.extend(valleys)
            valleys = []
            pid += 1
            values.append(PeakMatch(pid, p2, PeakFound.FOUND_PEAK, idx))
            p1 = p2
            scale = 1
            idx += 1
        elif dist < sed:
            # closer than expected, skip this peak
            idx += 1
        else:
            # further than expected, a valley before this peak
            pid += 1
            valleys.append(PeakMatch(pid, p1 + sed, PeakFound.FOUND_VALLEY, None))
            scale += 1
            if scale > max_scale_jump:
                _logger.warning("peak %s not found within expected distance from peak %s", pid, values[-1].count)
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
