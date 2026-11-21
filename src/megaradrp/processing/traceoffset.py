#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

"""Offset of the traces of an image relative to a TraceMap or a ModelMap

The traces of an image move with respect to the traces of the calibration
(flexures, temperature, time between exposures). This module measures that
displacement as the polynomial in the position of the trace along the slit
that ``extraction_offset`` takes: the offset of each trace is
``c0 + c1 * y``, with ``y`` the position of the trace in the reference column
of the calibration.

The image is collapsed in bands of columns. The profile of each band is
modelled as one Gaussian per trace, with a free amplitude, plus a smooth
background. For a given offset the model is linear, so its chi-square is
obtained by weighted least squares. The offset is the one with the lowest
chi-square summed over the bands. With free amplitudes the measurement does
not depend on which fibers are illuminated: the faint fibers have little
weight, but they do not bias it.

The distance between fibers is about 6 pixels, and the offsets can be of the
same order, so a periodic pattern alone would give aliases. The gaps between
the blocks of fibers and the ends of the slit break that degeneracy; the
contrast in chi-square between the best offset and the best alias, about one
fiber away, measures the reliability of the result.
"""

import dataclasses
import logging
import math
import warnings

import numpy
import scipy.optimize
import scipy.sparse
import scipy.sparse.linalg
import scipy.special
from scipy.interpolate import BSpline

_logger = logging.getLogger(__name__)

# a fiber with signal has an amplitude above this signal to noise ratio
SIGNAL_SNR = 5.0


@dataclasses.dataclass
class TraceOffset:
    """Offset of the traces of an image, and its diagnostics

    Attributes
    ----------
    coef : list of float
        [c0, c1], the coefficients of ``extraction_offset``
    error : list of float
        Uncertainty of `coef`, scaled with the reduced chi-square
    reliable : bool
        The offset is well defined: it is not at the edge of the search,
        there are fibers with signal, and the best alias is clearly worse
    alias_contrast : float
        Difference of chi-square between the best alias (more than half a
        fiber away) and the best offset, in units of the reduced chi-square
    alias_bootstrap : float
        Fraction of the bootstrap samples of the bands (drawn with
        replacement) where the best offset has a lower total chi-square than
        the best alias; with noise alone, about one half
    nsignal : int
        Number of fibers with signal
    sigma_scale : float
        Scale of the width of the profiles (the width in pixels for a TraceMap)
    reduced_chi2 : float
        Reduced chi-square of the best fit
    bands : list of dict
        For each band, its columns (x0, x1), its offset c0 (with c1 fixed) and
        the uncertainty
    grid, curve : numpy.ndarray
        The values of c0 of the global search and their chi-square (c1 = 0)
    message : str
        Why the result is not reliable, empty otherwise
    """

    coef: list
    error: list
    reliable: bool
    alias_contrast: float
    alias_bootstrap: float
    nsignal: int
    sigma_scale: float
    reduced_chi2: float
    bands: list
    grid: numpy.ndarray
    curve: numpy.ndarray
    message: str = ""


def _pixel_gaussian(rows, mu, sigma):
    """Gaussian of unit area integrated over pixels centred at `rows`"""
    scale = math.sqrt(2.0) * sigma
    upper = scipy.special.erf((rows + 0.5 - mu) / scale)
    lower = scipy.special.erf((rows - 0.5 - mu) / scale)
    return 0.5 * (upper - lower)


def _background_basis(rows, spacing):
    """Cubic B-spline basis on `rows`, with knots every `spacing` pixels"""
    lo, hi = rows[0] - 0.5, rows[-1] + 0.5
    nint = max(1, int(round((hi - lo) / spacing)))
    inner = numpy.linspace(lo, hi, nint + 1)
    knots = numpy.r_[[lo] * 3, inner, [hi] * 3]
    return BSpline.design_matrix(rows.astype(float), knots, 3)


class _Band:
    """Profile of a band of columns and its model"""

    def __init__(self, data, x0, x1, traces, ref_column, global_offset, base_sigma, ron, bg_spacing, mask=None):
        self.x0, self.x1 = x0, x1
        block = numpy.asarray(data[:, x0:x1], dtype=float)
        if mask is not None:
            block = numpy.where(mask[:, x0:x1], numpy.nan, block)
        ncols = numpy.count_nonzero(numpy.isfinite(block), axis=1)
        with warnings.catch_warnings():
            # rows with all their pixels masked
            warnings.simplefilter("ignore", RuntimeWarning)
            profile = numpy.nanmedian(block, axis=1)
        # position of each trace in the band: its mean over the columns, and its spread
        xs = numpy.linspace(x0, x1 - 1, 5)
        ys = numpy.array([[trace.polynomial(x) for x in xs] for trace in traces])
        self.y_ref = numpy.array([trace.polynomial(ref_column) for trace in traces])
        self.mu0 = ys.mean(axis=1) + global_offset(self.y_ref)
        spread = ys.std(axis=1)
        self.base_sigma = numpy.sqrt(numpy.asarray(base_sigma, dtype=float) ** 2 + spread**2)
        self.halfwidth = int(math.ceil(6 * self.base_sigma.max())) + 10
        ny = data.shape[0]
        lo = max(0, int(math.floor(self.mu0.min())) - 20)
        hi = min(ny, int(math.ceil(self.mu0.max())) + 21)
        self.rows = numpy.arange(lo, hi)
        good = numpy.isfinite(profile[lo:hi]) & (ncols[lo:hi] > 0)
        self.profile = numpy.where(good, profile[lo:hi], 0.0)
        # variance of the median of the columns
        var = 0.5 * math.pi * (ron**2 + numpy.clip(self.profile, 0, None)) / numpy.maximum(ncols[lo:hi], 1)
        self.weight = numpy.where(good, 1.0 / var, 0.0)
        self.background = _background_basis(self.rows, bg_spacing).tocsc()
        self.nfibers = len(traces)

    def design(self, offsets, sigma_scale):
        """Sparse design matrix: one Gaussian per trace, then the background"""
        mu = self.mu0 + offsets
        sigma = self.base_sigma * sigma_scale
        steps = numpy.arange(-self.halfwidth, self.halfwidth + 1)
        centre = numpy.rint(mu).astype(int)
        rows = centre[:, numpy.newaxis] + steps[numpy.newaxis, :]
        values = _pixel_gaussian(rows, mu[:, numpy.newaxis], sigma[:, numpy.newaxis])
        cols = numpy.repeat(numpy.arange(self.nfibers), steps.size).reshape(rows.shape)
        local = rows - self.rows[0]
        inside = (local >= 0) & (local < self.rows.size)
        gauss = scipy.sparse.csc_matrix(
            (values[inside], (local[inside], cols[inside])), shape=(self.rows.size, self.nfibers)
        )
        return scipy.sparse.hstack([gauss, self.background], format="csc")

    def fit(self, offsets, sigma_scale):
        """chi-square, amplitudes and design matrix for the offsets of the traces"""
        design = self.design(offsets, sigma_scale)
        weighted = design.multiply(self.weight[:, numpy.newaxis]).tocsc()
        normal = (design.T @ weighted).tocsc()
        rhs = weighted.T @ self.profile
        diag = normal.diagonal()
        ridge = 1e-10 * (diag.mean() if diag.size else 1.0)
        normal = normal + scipy.sparse.identity(normal.shape[0], format="csc") * ridge
        coeffs = scipy.sparse.linalg.spsolve(normal, rhs)
        resid = self.profile - design @ coeffs
        chi2 = float(numpy.sum(self.weight * resid**2))
        return chi2, coeffs, design

    def amplitudes(self, coeffs, design):
        """Amplitude of each trace and its inverse variance (neglecting the overlap of the profiles)"""
        gauss = design[:, : self.nfibers]
        ivar = numpy.asarray(gauss.multiply(gauss).T @ self.weight).ravel()
        return coeffs[: self.nfibers], ivar

    @property
    def npix(self):
        return int(numpy.count_nonzero(self.weight))


def valid_traces(trace_repr):
    """The traces of `trace_repr` with a valid position"""
    return [trace for trace in trace_repr.contents if trace.valid]


def _base_sigma(traces, columns):
    """Width of the profile of each trace in the columns, from a ModelMap; 1 for a TraceMap"""
    widths = []
    for trace in traces:
        model = getattr(trace, "model", None)
        stddev = model["params"].get("stddev") if model else None
        if stddev is None:
            widths.append(1.0)
        else:
            widths.append(float(numpy.mean([stddev(x) for x in columns])))
    return numpy.array(widths)


def _parabola_minimum(xs, ys):
    """Vertex and curvature of the parabola through the points around the minimum"""
    i = int(numpy.argmin(ys))
    lo, hi = max(0, i - 2), min(len(xs), i + 3)
    if hi - lo < 3:
        return xs[i], 0.0
    a, b, _ = numpy.polyfit(xs[lo:hi], ys[lo:hi], 2)
    if a <= 0:
        return xs[i], 0.0
    vertex = -b / (2 * a)
    return float(numpy.clip(vertex, xs[lo], xs[hi - 1])), float(a)


def measure_extraction_offset(
    data,
    trace_repr,
    nbands=16,
    band_width=64,
    search=(-8.0, 8.0),
    step=0.25,
    ron=3.4,
    mask=None,
    fit_slope=True,
    bg_spacing=200.0,
    alias_threshold=25.0,
    alias_bootstrap_threshold=0.99,
    min_signal=3,
):
    """Offset of the traces of `trace_repr` in the image `data`

    Parameters
    ----------
    data : numpy.ndarray
        Image after the basic processing (bias, dark, gain), in electrons
    trace_repr : TraceMap or ModelMap
        The traces of the calibration; its global_offset is applied
    nbands : int
        Number of bands of columns
    band_width : int
        Columns of each band
    search : tuple of float
        Range of the global search of c0, in pixels
    step : float
        Step of the global search, in pixels
    ron : float
        Readout noise, in electrons
    mask : numpy.ndarray or None
        Boolean array of the shape of `data`, True in the bad pixels
    fit_slope : bool
        Fit c1 when the fibers with signal cover enough of the slit
    bg_spacing : float
        Distance between the knots of the background, in pixels
    alias_threshold : float
        Minimum alias contrast of a reliable result
    alias_bootstrap_threshold : float
        Minimum fraction of the bootstrap samples of the bands where the best
        offset has a lower chi-square than the best alias
    min_signal : int
        Minimum number of fibers with signal of a reliable result

    Returns
    -------
    TraceOffset
    """
    traces = valid_traces(trace_repr)
    if not traces:
        raise ValueError("there are no valid traces")
    ncols = data.shape[1]
    start = numpy.percentile([max(t.start, 0) for t in traces], 90)
    stop = numpy.percentile([min(t.stop, ncols - 1) for t in traces], 10)
    centres = numpy.linspace(start + band_width, stop - band_width, nbands)
    bands = []
    for c in centres:
        x0 = int(round(c - band_width / 2))
        x1 = x0 + band_width
        base = _base_sigma(traces, numpy.linspace(x0, x1 - 1, 3))
        bands.append(
            _Band(data, x0, x1, traces, trace_repr.ref_column, trace_repr.global_offset, base, ron, bg_spacing, mask)
        )
    y_ref = bands[0].y_ref
    y_mid = 0.5 * (y_ref.min() + y_ref.max())
    is_model = any(getattr(t, "model", None) for t in traces)

    def total(c0, c1, scale, which=None):
        offsets = c0 + c1 * y_ref
        selected = bands if which is None else [bands[which]]
        return sum(band.fit(offsets, scale)[0] for band in selected)

    grid = numpy.arange(search[0], search[1] + 0.5 * step, step)

    def scan(scale):
        return numpy.array([total(c0, 0.0, scale) for c0 in grid])

    # width of the profiles: a scale of the ModelMap widths, or the width of a TraceMap in pixels
    scale = 1.0 if is_model else 1.5
    curve = scan(scale)
    c0 = grid[int(numpy.argmin(curve))]
    scales = numpy.arange(0.7, 1.61, 0.1) if is_model else numpy.arange(0.8, 3.01, 0.2)
    chi_scales = [total(c0, 0.0, s) for s in scales]
    best_scale = float(scales[int(numpy.argmin(chi_scales))])
    if abs(best_scale - scale) > 1e-6:
        scale = best_scale
        curve = scan(scale)
    ibest = int(numpy.argmin(curve))
    c0 = grid[ibest]

    # refine c0 around the minimum of the search
    fine = numpy.linspace(c0 - 2 * step, c0 + 2 * step, 21)
    chi_fine = numpy.array([total(x, 0.0, scale) for x in fine])
    c0, _ = _parabola_minimum(fine, chi_fine)
    c1 = 0.0

    # fibers with signal, and the part of the slit that they cover
    # the amplitudes of the bands combined with their weights: without signal, a normal distribution
    offsets = c0 + c1 * y_ref
    weighted_sum = numpy.zeros(len(traces))
    ivar_sum = numpy.zeros(len(traces))
    for band in bands:
        _, coeffs, design = band.fit(offsets, scale)
        amp, ivar = band.amplitudes(coeffs, design)
        weighted_sum += amp * ivar
        ivar_sum += ivar
    snr = weighted_sum / numpy.sqrt(numpy.maximum(ivar_sum, 1e-300))
    signal = snr > SIGNAL_SNR
    nsignal = int(numpy.count_nonzero(signal))
    span_all = numpy.ptp(y_ref)
    # the slope needs fibers with signal in both ends of the slit
    quarters = numpy.digitize(y_ref[signal], numpy.percentile(y_ref, [25, 50, 75]))
    per_quarter = numpy.bincount(quarters, minlength=4)
    slope_ok = per_quarter[0] >= min_signal and per_quarter[3] >= min_signal

    if fit_slope and slope_ok:
        # c1 as the difference of offset between the ends of the slit, to scale both parameters alike
        half = 0.5 * span_all

        def objective(p):
            return total(p[0] - p[1] * y_mid / half, p[1] / half, scale)

        start_point = numpy.array([c0 + c1 * y_mid, c1 * half])
        simplex = numpy.array([start_point, start_point + [0.3, 0.0], start_point + [0.0, 0.3]])
        res = scipy.optimize.minimize(
            objective,
            start_point,
            method="Nelder-Mead",
            options={"initial_simplex": simplex, "xatol": 0.005, "fatol": 0.01},
        )
        mid, delta = res.x
        c1 = delta / half
        c0 = mid - c1 * y_mid
    else:
        mid = c0

    chi_min = total(c0, c1, scale)
    npix = sum(band.npix for band in bands)
    nparams = sum(band.nfibers + band.background.shape[1] for band in bands) + 2
    reduced = chi_min / max(npix - nparams, 1)

    # uncertainty from the curvature of the chi-square, scaled with the reduced chi-square
    h0 = 0.05
    cp, cm = total(c0 + h0, c1, scale), total(c0 - h0, c1, scale)
    curv0 = (cp + cm - 2 * chi_min) / h0**2
    err0 = math.sqrt(2 * reduced / curv0) if curv0 > 0 else math.inf
    err1 = 0.0
    if c1 != 0.0:
        h1 = h0 / (0.5 * span_all)
        cp1, cm1 = total(c0 - h1 * y_mid, c1 + h1, scale), total(c0 + h1 * y_mid, c1 - h1, scale)
        curv1 = (cp1 + cm1 - 2 * chi_min) / h1**2
        err1 = math.sqrt(2 * reduced / curv1) if curv1 > 0 else math.inf

    # alias contrast: the best value of the search more than half a fiber away from its minimum
    gaps = numpy.diff(numpy.sort(y_ref))
    pitch = float(numpy.median(gaps[gaps < 3 * numpy.median(gaps)]))
    far = numpy.abs(grid - grid[ibest]) >= 0.5 * pitch
    if far.any():
        ialias = numpy.flatnonzero(far)[int(numpy.argmin(curve[far]))]
        alias = (curve[ialias] - curve[ibest]) / reduced
        # the difference must be consistent between the bands, not the result of a few of them
        diffs = numpy.array(
            [
                total(grid[ialias], 0.0, scale, which=i) - total(grid[ibest], 0.0, scale, which=i)
                for i in range(len(bands))
            ]
        )
        # the total must favour the best offset whichever bands are used: bootstrap over the bands
        rng = numpy.random.default_rng(0)
        resampled = diffs[rng.integers(0, len(diffs), size=(2000, len(diffs)))].sum(axis=1)
        alias_bootstrap = float(numpy.mean(resampled > 0))
    else:
        alias, alias_bootstrap = math.inf, 1.0

    # offset of each band, with c1 fixed
    band_results = []
    for i, band in enumerate(bands):
        local = numpy.linspace(mid - 1.0, mid + 1.0, 21)
        chi_local = numpy.array([total(x - c1 * y_mid, c1, scale, which=i) for x in local])
        vertex, a = _parabola_minimum(local, chi_local)
        band_reduced = chi_local.min() / max(band.npix - band.nfibers - band.background.shape[1], 1)
        band_err = math.sqrt(band_reduced / a) if a > 0 else math.inf
        band_results.append({"x0": band.x0, "x1": band.x1, "c0": vertex - c1 * y_mid, "error": band_err})

    messages = []
    if ibest in (0, len(grid) - 1):
        messages.append("the minimum is at the edge of the search")
    if nsignal < min_signal:
        messages.append(f"only {nsignal} fibers with signal")
    if not alias > alias_threshold:
        messages.append(f"alias contrast {alias:.1f} below {alias_threshold}")
    if alias_bootstrap < alias_bootstrap_threshold:
        messages.append(f"the offset is preferred to the alias in {alias_bootstrap:.1%} of the bootstrap samples")
    result = TraceOffset(
        coef=[float(c0), float(c1)],
        error=[err0, err1],
        reliable=not messages,
        alias_contrast=float(alias),
        alias_bootstrap=alias_bootstrap,
        nsignal=nsignal,
        sigma_scale=scale,
        reduced_chi2=float(reduced),
        bands=band_results,
        grid=grid,
        curve=curve,
        message="; ".join(messages),
    )
    _logger.info(
        "extraction offset [%.4f, %.4g] (+- %.4f, %.2g), alias contrast %.1f, %d fibers with signal%s",
        result.coef[0],
        result.coef[1],
        err0,
        err1,
        alias,
        nsignal,
        "" if result.reliable else ", not reliable: " + result.message,
    )
    return result


def combine_extraction_offsets(results, trace_repr, max_spread=0.3, logger=None):
    """Offset of several images of an observation, from the measurement in each one

    The offset is the mean of the coefficients of the reliable measurements.
    If the offsets of the reliable measurements in the middle of the slit
    differ by more than `max_spread` pixels, there is a warning: the images
    are combined before the extraction, so a single offset is applied to all
    of them. Without reliable measurements the offset is 0, with a warning.

    Parameters
    ----------
    results : list of TraceOffset
        The measurement in each image
    trace_repr : TraceMap or ModelMap
        The traces of the calibration
    max_spread : float
        Maximum difference, in pixels, between the offsets of the images
    logger : logging.Logger, optional

    Returns
    -------
    coef : list of float
        The coefficients of ``extraction_offset``
    history : list of str
        Lines for the HISTORY of the product
    """
    logger = logger or _logger
    y_ref = numpy.array([trace.polynomial(trace_repr.ref_column) for trace in valid_traces(trace_repr)])
    y_mid = 0.5 * (y_ref.min() + y_ref.max())
    history = []
    for idx, result in enumerate(results, start=1):
        status = "reliable" if result.reliable else f"not reliable ({result.message})"
        history.append(f"Extraction offset of image {idx}: [{result.coef[0]:.4f}, {result.coef[1]:.4g}], {status}")
    reliable = [result for result in results if result.reliable]
    if not reliable:
        logger.warning("no reliable measurement of the extraction offset in %d images, using 0", len(results))
        history.append("No reliable measurement of the extraction offset, using 0")
        return [0.0], history
    coef = numpy.mean([result.coef for result in reliable], axis=0).tolist()
    middle = [result.coef[0] + result.coef[1] * y_mid for result in reliable]
    spread = max(middle) - min(middle)
    if spread > max_spread:
        logger.warning(
            "the extraction offsets of the images differ by %.2f pixels in the middle of the slit, "
            "more than %.2f; their mean is applied to the combined image",
            spread,
            max_spread,
        )
        history.append(f"The extraction offsets of the images differ by {spread:.2f} pixels")
    if len(reliable) < len(results):
        logger.warning(
            "the extraction offset is reliable in %d of %d images, using their mean", len(reliable), len(results)
        )
    logger.info("extraction offset [%.4f, %.4g] from %d images", coef[0], coef[1], len(reliable))
    history.append(f"Extraction offset measured: [{coef[0]:.4f}, {coef[1]:.4g}]")
    return coef, history
