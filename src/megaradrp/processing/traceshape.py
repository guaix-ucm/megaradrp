#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

"""Shape of the fiber traces predicted from the neighbouring fibers

The traces of neighbouring fibers have almost the same shape, and the shape
changes slowly along the pseudoslit. The shape of the complete traces, as a
function of their position in the reference column, predicts the shape of
any trace. It is used to correct the traced points that jump to a
neighbouring fiber and to complete the traces that stop before the borders
of the image, where the fibers are faint or defocused.
"""

import logging

import numpy
import numpy.polynomial.polynomial as nppol

from megaradrp.products.tracemap import GeometricTrace

_logger = logging.getLogger(__name__)


class TraceShapeModel:
    """Shape of the traces, y(x) - y(xref), as a function of y(xref)

    For each column of `xgrid`, the shapes are fitted with a polynomial of
    degree `deg` in the position of the trace in the reference column, with
    an iterative rejection of outliers.

    Parameters
    ----------
    xgrid : array_like
        Columns where the shapes are sampled
    yref : array_like
        Positions of the traces in the reference column
    shapes : array_like
        Shapes of the traces, y(xgrid) - yref, with shape (len(yref), len(xgrid))
    deg : int
        Degree of the polynomial in yref
    nsig : float
        Rejection threshold, in units of the robust standard deviation
    niter : int
        Maximum number of rejection iterations
    min_sigma : float
        Residuals below this value (pixels) are never rejected
    """

    def __init__(self, xgrid, yref, shapes, deg=4, nsig=4.0, niter=5, min_sigma=0.05):
        self.xgrid = numpy.asarray(xgrid, dtype=float)
        yref = numpy.asarray(yref, dtype=float)
        shapes = numpy.asarray(shapes, dtype=float)
        if shapes.shape != (len(yref), len(self.xgrid)):
            raise ValueError("shapes must have shape (len(yref), len(xgrid))")
        if len(yref) <= deg:
            raise ValueError(f"at least {deg + 1} traces are required to fit a polynomial of degree {deg}")
        self.center = 0.5 * (yref.min() + yref.max())
        self.scale = max(0.5 * (yref.max() - yref.min()), 1.0)
        u = (yref - self.center) / self.scale
        coefs = numpy.empty((deg + 1, len(self.xgrid)))
        for j in range(len(self.xgrid)):
            keep = numpy.ones(len(u), dtype=bool)
            for _ in range(niter):
                c = nppol.polyfit(u[keep], shapes[keep, j], deg)
                res = shapes[:, j] - nppol.polyval(u, c)
                sigma = 1.4826 * numpy.median(numpy.abs(res[keep]))
                new_keep = numpy.abs(res) <= max(nsig * sigma, min_sigma)
                if numpy.array_equal(new_keep, keep) or new_keep.sum() <= deg:
                    break
                keep = new_keep
            coefs[:, j] = c
        self.coefs = coefs

    def shape(self, yref, x):
        """Shape y(x) - yref of the trace with position yref in the reference column"""
        u = (yref - self.center) / self.scale
        values = nppol.polyval(u, self.coefs)
        return numpy.interp(x, self.xgrid, values)


def _crossing_traces(ids, ys):
    """Fibers whose trace comes closer than half the usual distance to the next fiber"""
    crossing = set()
    for a, b in zip(ids[:-1], ids[1:]):
        if b == a + 1:
            sep = ys[b] - ys[a]
            if sep.min() < 0.5 * numpy.median(sep):
                crossing.update((a, b))
    return crossing


def complete_traces(
    traces, points, xlimits, xref, poldeg, step=20, deg=4, max_dev=1.5, ref_window=300, min_complete=20
):
    """Correct and complete the traces with the shape of the neighbouring traces

    The shape model is fitted with the complete traces: traced from `xlimits[0]`
    to `xlimits[1]` and not crossing a neighbour. For each trace, the offset to
    the model is measured with the points within `ref_window` pixels of the
    reference column, where the tracing starts. The points farther than
    `max_dev` pixels from the model plus the offset are rejected. A trace that
    is not complete or has rejected points is fitted again with its accepted
    points, and with the model plus the offset beyond them, over the whole
    range.

    Parameters
    ----------
    traces : list of GeometricTrace
        The traces, changed in place
    points : dict
        Traced points of each fiber, fibid: (x, y)
    xlimits : tuple of int
        First and last columns of a complete trace
    xref : float
        Reference column, where the tracing starts
    poldeg : int
        Degree of the polynomial of the traces
    step : int
        Distance between the columns where the shapes are sampled
    deg : int
        Degree of the polynomial of the shape model in the position of the traces
    max_dev : float
        Maximum distance (pixels) of a traced point to the model
    ref_window : float
        Half width (pixels) of the region around `xref` used to measure the offset
    min_complete : int
        Minimum number of complete traces required to fit the model

    Returns
    -------
    list of int
        The fibids of the traces fitted again
    """
    xstart, xend = xlimits
    xgrid = numpy.arange(xstart, xend + step, step, dtype=float)
    xgrid[-1] = min(xgrid[-1], xend)

    index = {t.fibid: idx for idx, t in enumerate(traces) if t.valid and t.fibid in points}
    ids = sorted(index)
    if not ids:
        return []
    ys = {f: traces[index[f]].polynomial(xgrid) for f in ids}
    crossing = _crossing_traces(ids, ys)

    def is_complete(fibid):
        x = points[fibid][0]
        return x.min() <= xstart + step and x.max() >= xend - step

    complete = [f for f in ids if is_complete(f) and f not in crossing]
    if len(complete) < min_complete:
        _logger.warning("only %d complete traces, the traces are not completed with a model", len(complete))
        return []

    yref = numpy.array([traces[index[f]].polynomial(xref) for f in complete])
    shapes = numpy.array([ys[f] - traces[index[f]].polynomial(xref) for f in complete])
    model = TraceShapeModel(xgrid, yref, shapes, deg=deg)
    _logger.debug("shape model fitted with %d complete traces", len(complete))

    changed = []
    for fibid in ids:
        trace = traces[index[fibid]]
        x, y = points[fibid]
        # the position in the reference column, from the traced points
        y0 = y[numpy.argmin(numpy.abs(x - xref))]
        res = y - (y0 + model.shape(y0, x))
        near = numpy.abs(x - xref) <= ref_window
        offset = numpy.median(res[near]) if near.sum() >= 5 else numpy.median(res)
        good = numpy.abs(res - offset) <= max_dev
        if is_complete(fibid) and good.all():
            continue
        if good.sum() < poldeg + 1:
            _logger.warning("fiber %d has only %d points close to the model, trace not changed", fibid, good.sum())
            continue
        xg = x[good]
        yg = y[good]
        outside = (xgrid < xg.min()) | (xgrid > xg.max())
        xfit = numpy.concatenate([xg, xgrid[outside]])
        yfit = numpy.concatenate([yg, y0 + offset + model.shape(y0, xgrid[outside])])
        pfit = nppol.polyfit(xfit, yfit, poldeg)
        traces[index[fibid]] = GeometricTrace(
            fibid=trace.fibid, boxid=trace.boxid, start=int(xstart), stop=int(xend), fitparms=pfit.tolist()
        )
        changed.append(fibid)
        _logger.debug(
            "fiber %d: traced %d-%d, %d of %d points rejected, completed with the shape model",
            fibid,
            x.min(),
            x.max(),
            (~good).sum(),
            len(good),
        )
    return changed
