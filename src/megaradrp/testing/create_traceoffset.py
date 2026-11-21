"""Synthetic images with the traces moved, to test the measurement of extraction_offset"""

import numpy

from megaradrp.processing.traceoffset import _pixel_gaussian
from megaradrp.products.tracemap import GeometricTrace, TraceMap

NROWS, NCOLS = 320, 640
PITCH = 6.0
SIGMA = 1.3


def fiber_positions():
    """Position of the fibers in the reference column: 4 blocks of 10 fibers, with gaps"""
    positions = []
    for block in range(4):
        start = 30 + block * (10 * PITCH + 15)
        positions.extend(start + PITCH * numpy.arange(10) + 0.3 * numpy.sin(numpy.arange(10)))
    return numpy.array(positions)


def create_offset_tracemap():
    """TraceMap of 40 slightly tilted traces in 4 blocks, plus a missing fiber"""
    tracemap = TraceMap(instrument="MEGARA")
    tracemap.ref_column = NCOLS // 2
    tracemap.total_fibers = 41
    slope = 0.004
    for i, y in enumerate(fiber_positions()):
        fitparms = [y - slope * NCOLS / 2, slope]
        tracemap.contents.append(GeometricTrace(i + 1, i // 10 + 1, 0, NCOLS - 1, fitparms=fitparms))
    tracemap.contents.append(GeometricTrace(41, 4, 0, NCOLS - 1, fitparms=[]))
    return tracemap


def create_offset_image(tracemap, c0, c1=0.0, flux=None, sky=0.0, seed=0):
    """Image with the traces of `tracemap` moved by c0 + c1 * y, in electrons, with Poisson and readout noise

    `y` is the position of each trace in the reference column of `tracemap`.
    """
    rng = numpy.random.default_rng(seed)
    traces = [t for t in tracemap.contents if t.valid]
    flux = numpy.full(len(traces), 2000.0) if flux is None else numpy.asarray(flux, dtype=float)
    rows = numpy.arange(NROWS)[:, numpy.newaxis]
    cols = numpy.arange(NCOLS)
    image = numpy.full((NROWS, NCOLS), float(sky))
    for trace, f in zip(traces, flux):
        y_ref = trace.polynomial(tracemap.ref_column)
        mu = trace.polynomial(cols) + c0 + c1 * y_ref
        image += f * _pixel_gaussian(rows, mu[numpy.newaxis, :], SIGMA)
    return rng.poisson(image).astype(float) + rng.normal(0, 3.4, image.shape)
