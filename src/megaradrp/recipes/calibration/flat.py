#
# Copyright 2011-2025 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

"""Fiber flat calibration Recipe for Megara"""

import numpy
from astropy.io import fits
from numina.core import Result, Parameter
import numina.exceptions

from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.instrument.focalplane import FocalPlaneConf
from megaradrp.ntypes import MasterFiberFlat
import megaradrp.requirements as reqs
from megaradrp.ntypes import ProcessedRSS, ProcessedFrame

# Flat 2D
from numina.processing.combine import basic_processing_with_combination
from numina.array import combine

# Create RSS
from megaradrp.processing.aperture import ApertureExtractor
from megaradrp.processing.wavecalibration import WavelengthCalibrator
from megaradrp.processing.fiberflat import Splitter, FlipLR


def absorption_windows(resolution):
    """Windows of the absorptions in the light of the calibration lamps (air, Angstrom)

    The windows are computed for each resolution class ("LR", "MR", "HR") by
    tools/flat_absorption_windows.py: O2 and H2O lines along the path of the
    calibration light, and the Na D lines of the lamp.
    """
    import pkgutil

    from astropy.table import Table

    data = pkgutil.get_data("megaradrp.instrument.configs", "flat_absorption_windows.ecsv")
    table = Table.read(data.decode("utf8"), format="ascii.ecsv")
    sel = table["resolution"] == resolution
    return [(float(a), float(b)) for a, b in zip(table["wl_min"][sel], table["wl_max"][sel])]


def _flat_knots(x_fit, step, min_points=10):
    """Interior knots of the spline of the average spectrum

    Equally spaced between the ends of `x_fit`, not farther than `step`
    pixels. A knot with less than `min_points` pixels of `x_fit` since the
    previous one is removed.
    """
    first, last = x_fit[0], x_fit[-1]
    nknots = max(int(numpy.ceil((last - first) / step)), 1)
    knots = numpy.linspace(first, last, nknots + 1)[1:-1]
    kept, prev = [], first
    for knot in knots:
        if numpy.count_nonzero((x_fit > prev) & (x_fit <= knot)) >= min_points:
            kept.append(knot)
            prev = knot
    while kept and numpy.count_nonzero(x_fit > kept[-1]) < min_points:
        kept.pop()
    return kept


def _extend_linear(spline, x, first, last):
    """Evaluate `spline` in `x`, continued beyond [first, last] with straight lines

    Each line has the value and the slope of the spline at the border.
    """
    out = spline(x)
    deriv = spline.derivative()
    for border, side in ((first, x < first), (last, x > last)):
        if numpy.any(side):
            out[side] = spline(border) + deriv(border) * (x[side] - border)
    return out


def interpolate_windows(data, skip, coverage, nside=5):
    """Replace the pixels in `skip` with a straight line between both sides, in each row

    Each side is the median of the `nside` closest pixels outside `skip` with
    coverage. If one side is missing, its value is the other side.

    Parameters
    ----------
    data : numpy.ndarray
        Array of shape (nrows, npix), changed in place
    skip : numpy.ndarray
        Boolean array of shape (npix,), the pixels to replace
    coverage : numpy.ndarray
        Boolean array of shape (nrows, npix), the pixels with data
    """
    npix = data.shape[1]
    edges = numpy.flatnonzero(numpy.diff(numpy.r_[0, skip.astype(int), 0]))
    for a, b in zip(edges[::2], edges[1::2]):
        left = numpy.arange(max(0, a - 3 * nside), a)
        right = numpy.arange(b, min(npix, b + 3 * nside))
        left = left[~skip[left]][-nside:]
        right = right[~skip[right]][:nside]
        sides = []
        for idx in (left, right):
            if idx.size:
                with numpy.errstate(all="ignore"):
                    sides.append(numpy.nanmedian(numpy.where(coverage[:, idx], data[:, idx], numpy.nan), axis=1))
            else:
                sides.append(None)
        vleft, vright = sides
        if vleft is None and vright is None:
            continue
        if vleft is None:
            vleft = vright
        if vright is None:
            vright = vleft
        vleft = numpy.where(numpy.isfinite(vleft), vleft, vright)
        vright = numpy.where(numpy.isfinite(vright), vright, vleft)
        # the sides are placed at the centre of their pixels
        xl = numpy.median(left) if left.size else a - 1
        xr = numpy.median(right) if right.size else b
        frac = (numpy.arange(a, b) - xl) / max(xr - xl, 1)
        line = vleft[:, numpy.newaxis] + (vright - vleft)[:, numpy.newaxis] * frac[numpy.newaxis, :]
        ok = numpy.isfinite(line)
        block = data[:, a:b]
        block[ok] = line[ok]


# Flags of a fiber in the fiber flat that make it invalid
INVALID_FLAT_FLAGS = frozenset(["broken"])


def _trace_polynomials(tracemap):
    """(fibid, position) of the valid traces of a TraceMap or a ModelMap, with its global offset

    The position is a function giving the row of the centre of the fiber for
    each column: the polynomial of a TraceMap, or the spline of the centre
    of the profiles of a ModelMap.
    """
    out = []
    for trace in tracemap.contents:
        if not trace.valid:
            continue
        try:
            centre = trace.polynomial
        except (KeyError, TypeError):
            # a fiber of a ModelMap without model
            continue
        offset = float(tracemap.global_offset(centre(tracemap.ref_column)))

        def position(x, centre=centre, offset=offset):
            return centre(x) + offset

        out.append((trace.fibid, position))
    return out


def find_zero_order(image, traces, threshold=0.02, min_fwhm=5.0, max_fwhm=40.0, width=151, nbins=10):
    """Find the zero order of the VPH in a 2D flat image

    The zero order is diffuse light added at a fixed range of columns, also
    between the fibers. In each row midway between neighbouring traces the
    image is divided by the level of the fibers next to it; the median over
    groups of those rows, minus its running median along the columns
    (`width` pixels), keeps the narrow excess. A band is accepted if its
    peak is above `threshold` (fraction of the level of the fibers) and its
    FWHM is between `min_fwhm` and `max_fwhm` pixels. Weaker bands (1-2 %
    in some gratings) are left: interpolating across them would also remove
    the real structure of the flat there, such as the etalon of the BLUE
    order sorting filter, by about as much.

    The band is slightly tilted. Its centre is a straight line in the row,
    fitted to the centres of Gaussians fitted in each group of rows within
    +-3 FWHM of the peak; groups where the fit fails, the band is weaker
    than half the threshold or its FWHM is far from that of the band, and
    outliers of the line (3 sigma), are left out.

    Parameters
    ----------
    image : numpy.ndarray
        2D image of the flat, before the extraction
    traces : list of (int, callable)
        Fiber id and position (row) of each valid trace as a function of the column

    Returns
    -------
    list of dict
        For each band: `center`, a numpy Polynomial of degree 1 giving its
        column as a function of the row; `fwhm` in pixels; `amplitude`, the peak excess
        between the fibers relative to the level of the fibers; `rows`, the
        range of rows where it was measured.
    """
    from scipy.ndimage import median_filter
    from scipy.optimize import curve_fit
    from scipy.signal import find_peaks, peak_widths

    def gauss(x, amp, centre, sigma, base):
        return amp * numpy.exp(-0.5 * ((x - centre) / sigma) ** 2) + base

    nrows, ncols = image.shape
    xx = numpy.arange(ncols)
    if len(traces) < 2:
        return []
    ypos = numpy.array([poly(xx) for _, poly in traces])
    order = numpy.argsort(ypos[:, ncols // 2])
    ypos = ypos[order]
    gaps = numpy.diff(ypos[:, ncols // 2])
    near = gaps < 1.5 * numpy.median(gaps)
    mid = numpy.rint(0.5 * (ypos[:-1] + ypos[1:])[near]).astype(int)
    lo_rows = numpy.rint(ypos[:-1][near]).astype(int)
    hi_rows = numpy.rint(ypos[1:][near]).astype(int)
    for rows in (mid, lo_rows, hi_rows):
        numpy.clip(rows, 0, nrows - 1, out=rows)
    valleys = image[mid, xx]
    centre = slice(ncols // 2 - 100, ncols // 2 + 100)
    level = 0.5 * (
        numpy.median(image[lo_rows, xx][:, centre], axis=1) + numpy.median(image[hi_rows, xx][:, centre], axis=1)
    )
    use = level > 0.3 * numpy.median(level)
    if numpy.count_nonzero(use) < nbins:
        return []
    valleys = valleys[use] / level[use, numpy.newaxis]
    row_of = mid[use, ncols // 2]
    groups = numpy.array_split(numpy.arange(valleys.shape[0]), nbins)
    excess = []
    for group in groups:
        prof = numpy.median(valleys[group], axis=0)
        excess.append(prof - median_filter(prof, size=width, mode="nearest"))
    excess = numpy.array(excess)
    total = numpy.median(excess, axis=0)
    total[: width // 2] = 0
    total[-(width // 2) :] = 0
    peaks, _ = find_peaks(total, height=threshold, width=(min_fwhm, max_fwhm))
    bands = []
    for peak in peaks:
        fwhm = float(peak_widths(total, [peak], rel_height=0.5)[0][0])
        half = int(numpy.ceil(3 * fwhm))
        window = numpy.arange(max(peak - half, 0), min(peak + half + 1, ncols))
        sigma0 = fwhm / 2.3548
        rows, cols = [], []
        for group, prof in zip(groups, excess):
            seg = prof[window]
            try:
                (amp, centre, sigma, _), _ = curve_fit(
                    gauss, window, seg, p0=(max(seg.max(), threshold), window[numpy.argmax(seg)], sigma0, 0.0)
                )
            except (RuntimeError, ValueError):
                continue
            fwhm_group = 2.3548 * abs(sigma)
            good = amp > 0.5 * threshold and min_fwhm <= fwhm_group <= max_fwhm and 0.5 < fwhm_group / fwhm < 1.6
            if good and window[0] < centre < window[-1]:
                rows.append(numpy.median(row_of[group]))
                cols.append(centre)
        rows, cols = numpy.array(rows), numpy.array(cols)
        keep = numpy.ones(rows.size, dtype=bool)
        if rows.size >= 3:
            for _ in range(3):
                line = numpy.polynomial.Polynomial.fit(rows[keep], cols[keep], 1)
                resid = cols - line(rows)
                scatter = max(1.4826 * numpy.median(numpy.abs(resid[keep])), 0.5)
                new_keep = numpy.abs(resid) < 3 * scatter
                if new_keep.sum() < 3 or numpy.array_equal(new_keep, keep):
                    break
                keep = new_keep
            center = numpy.polynomial.Polynomial.fit(rows[keep], cols[keep], 1)
        elif rows.size:
            center = numpy.polynomial.Polynomial([numpy.median(cols)])
            keep = numpy.ones(rows.size, dtype=bool)
        else:
            center = numpy.polynomial.Polynomial([float(peak)])
            rows = numpy.array([numpy.median(row_of)])
            keep = numpy.ones(1, dtype=bool)
        rows = rows[keep]
        bands.append(
            {
                "center": center,
                "fwhm": fwhm,
                "amplitude": float(total[peak]),
                "rows": (float(min(rows)), float(max(rows))),
            }
        )
    return bands


def remove_zero_order(rss, traces, bands, scale=1.5, nside=5):
    """Interpolate each fiber of an extracted RSS across the columns of the zero order

    `rss` is the extracted RSS before the flip and the wavelength calibration
    (row fibid - 1, detector columns), changed in place. In each fiber the
    columns within `scale` times the FWHM of the centre of each band, at the
    row of the fiber, are replaced by a straight line between the medians of
    the `nside` pixels at each side.
    """
    ncols = rss.shape[1]
    for band in bands:
        half = scale * band["fwhm"]
        col0 = band["center"](numpy.mean(band["rows"]))
        for fibid, poly in traces:
            col = band["center"](poly(col0))
            a = max(int(numpy.floor(col - half)), 0)
            b = min(int(numpy.ceil(col + half)) + 1, ncols)
            if a <= 0 or b >= ncols:
                continue
            row = rss[fibid - 1]
            left = numpy.median(row[max(a - nside, 0) : a])
            right = numpy.median(row[b : b + nside])
            xl = a - (nside + 1) / 2
            xr = b + (nside - 1) / 2
            row[a:b] = left + (right - left) * (numpy.arange(a, b) - xl) / (xr - xl)


def flag_fibers(flat, coverage, faint, good, col1, col2, skip=None, max_dev=0.3, min_frac=0.02, nneigh=5, width=51):
    """Quality flags of the fibers in a fiber flat

    - "faint": the fiber is fainter than the limit of the average spectrum
      (informational)
    - "deviant": the shape of the flat of the fiber differs by more than
      `max_dev` from the median shape of its `nneigh` neighbors at each side
      on the pseudo-slit, in more than `min_frac` of its pixels with coverage,
      after a median filter of `width` pixels (informational: in the
      commissioning data these are steep features of the spectrum with small
      differences in the wavelength calibration, and traces off-center where
      the signal is weak; both are also in the science spectra)
    - "broken": the level of the fiber is not positive, or the flat has non
      finite values where the fiber has coverage (the fiber is invalid)

    The shape of each fiber is its flat divided by its median in the columns
    `col1:col2`.

    Parameters
    ----------
    flat : numpy.ndarray
        The fiber flat, shape (nfibers, npix)
    coverage : numpy.ndarray
        Boolean array, the pixels with coverage
    faint : numpy.ndarray
        Boolean array of the fibers fainter than the limit
    good : numpy.ndarray
        Boolean array of the fibers not marked invalid before
    skip : numpy.ndarray, optional
        Boolean array of the pixels interpolated in the absorption windows

    Returns
    -------
    dict
        fibid: list of flags, for the fibers with flags
    """
    import warnings

    from scipy.ndimage import median_filter

    nfibers = flat.shape[0]
    use = coverage.copy()
    if skip is not None:
        use &= ~skip[numpy.newaxis, :]
    checked = good & coverage.any(axis=1)
    level = numpy.median(flat[:, col1:col2], axis=1)
    finite = numpy.all(numpy.isfinite(flat) | ~coverage, axis=1)
    broken = checked & ~((level > 0) & finite)
    shaped = checked & ~broken
    shape = numpy.full(flat.shape, numpy.nan)
    shape[shaped] = numpy.where(coverage[shaped], flat[shaped] / level[shaped, numpy.newaxis], numpy.nan)

    flags = {}
    for idx in numpy.flatnonzero(checked):
        fiber_flags = []
        if faint[idx]:
            fiber_flags.append("faint")
        if broken[idx]:
            fiber_flags.append("broken")
        else:
            neigh = [j for j in range(max(0, idx - nneigh), min(nfibers, idx + nneigh + 1)) if j != idx and shaped[j]]
            if neigh:
                with warnings.catch_warnings():
                    # pixels without coverage in all the neighbors
                    warnings.simplefilter("ignore", RuntimeWarning)
                    ratio = shape[idx] / numpy.nanmedian(shape[neigh], axis=0)
                pixels = use[idx] & numpy.isfinite(ratio)
                ratio = median_filter(numpy.where(numpy.isfinite(ratio), ratio, 1.0), width)
                if pixels.any() and numpy.mean(numpy.abs(ratio[pixels] - 1) > max_dev) > min_frac:
                    fiber_flags.append("deviant")
        if fiber_flags:
            flags[int(idx) + 1] = fiber_flags
    return flags


def check_qc_fiberflat(fibers_hdr, nfibers=None):
    """Global quality of the fiber flat

    GOOD if no fiber is invalid by the checks of the fiber flat, PARTIAL if some
    are, BAD if more than 5 % of the fibers are.
    """
    from numina.types.qc import QC

    if nfibers is None:
        nfibers = fibers_hdr.get("NFIBERS", 623)
    invalid = 0
    for key, value in fibers_hdr.items():
        if key.endswith("QF") and not INVALID_FLAT_FLAGS.isdisjoint(str(value).split(",")):
            invalid += 1
    if invalid > 0.05 * nfibers:
        return QC.BAD
    if invalid > 0:
        return QC.PARTIAL
    return QC.GOOD


def _smoothing_window_check(val):
    """Check 'smoothing_window'"""
    # raise ValueError("something")
    if val < 5:
        raise numina.exceptions.ValidationError("must be >= 5")
    if val % 2 == 0:
        raise numina.exceptions.ValidationError("must be odd")
    return val


class FiberFlatRecipe(MegaraBaseRecipe):
    """Process FIBER_FLAT images and create MASTER_FIBER_FLAT product.

    This recipe process a set of continuum flat images obtained in
    **Fiber Flat** mode and returns the master fiber flat product
    The recipe also returns the result of processing the input images up to
    slitflat correction. and the result RSS of the processing
    up to wavelength calibration.

    See Also
    --------
    megaradrp.products.MasterFiberFlat: description of MasterFiberFlat product
    megaradrp.processing.aperture: aperture extraction
    megaradrp.processing.wavecalibration: resampling for wavelength calibration

    Notes
    -----
    Images provided in `obresult` are trimmed and corrected from overscan,
    bad pixel mask (if `master_bpm` is not None), bias and dark current
    (if `master_dark` is not None) and corrected from pixel-to-pixel flat
    if `master_slitflat` is not None.
    Images thus corrected are the stacked using the median.

    The result of the combination is saved as an intermediate result, named
    'reduced_image.fits'. This combined image is also returned in the field
    `reduced_image` of the recipe result.

    The apertures in the 2D image are extracted, using the information in
    `master_apertures` and resampled accoding to the wavelength calibration in
    `master_wlcalib`. The resulting RSS is saved as an intermediate
    result named 'reduced_rss.fits'. This RSS is also returned in the field
    `reduced_rss` of the recipe result.

    To normalize the `master_fiberflat`, each fiber is divided by a smoothed
    version (using a Savitzky-Golay filter) of the average of the valid fibers.
    Finally, all the pixels with information are fiiled with ones. This RSS
    image is returned in the field `master_fiberflat` of the recipe result.

    """

    # Requirements
    method = Parameter("median", description="Combination method", choices=["mean", "median", "mediancr", "sigmaclip"])
    method_kwargs = Parameter(dict(), description="Arguments for combination method", optional=True)
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_bpm = reqs.MasterBPMRequirement()
    master_slitflat = reqs.MasterSlitFlatRequirement()
    master_apertures = reqs.MasterAperturesRequirement(alias="master_traces")
    smoothing_window = Parameter(31, "Window for smoothing (must be odd)", validator=_smoothing_window_check)
    skip_absorption_windows = Parameter(
        True, "Skip the windows of the absorptions of the calibration light (air and Na D), interpolating across them"
    )
    smoothing_knot_step = Parameter(
        150.0, "Distance in Angstrom between the knots of the spline that smooths the average spectrum of the fibers"
    )
    remove_zero_order = Parameter(
        True, "Find the zero order of the VPH (diffuse light at fixed columns) and interpolate each fiber across it"
    )
    extraction_offset = Parameter([0.0], "Offset traces for extraction", accept_scalar=True)
    master_wlcalib = reqs.WavelengthCalibrationRequirement()

    # Results
    reduced_image = Result(ProcessedFrame)
    reduced_rss = Result(ProcessedRSS)
    master_fiberflat = Result(MasterFiberFlat)

    def process_flat2d(self, rinput):
        flow = self.init_filters(rinput, rinput.obresult.configuration)
        fmethod = getattr(combine, rinput.method)
        final_image = basic_processing_with_combination(
            rinput,
            flow,
            method=fmethod,
            method_kwargs=rinput.method_kwargs,
        )
        hdr = final_image[0].header
        self.set_base_headers(hdr)
        return final_image

    def obtain_fiber_flat(
        self,
        rss_wl,
        col1=1900,
        col2=2100,
        window=31,
        degree=3,
        knot_step=150.0,
        min_level=0.3,
        full_coverage=0.99,
        skip_windows=None,
    ):
        """Divide each fiber by the smoothed average spectrum of the fibers

        Each valid fiber is normalized by its median in the columns `col1:col2`.
        The fibers fainter than `min_level` times the median of the fibers are
        left out. In each wavelength, the average spectrum is the mean of the
        normalized fibers with coverage (WLMAP), divided by the number of those
        fibers.

        The average spectrum is smoothed with a cubic spline with knots every
        `knot_step` Angstrom, fitted only where all the fibers have coverage
        (at least `full_coverage` times the largest number of fibers); a knot
        with less than 10 fitted pixels since the previous one is removed.
        Near the borders only some fibers have coverage, and their average
        changes with the fibers that enter it, so the spline is continued
        there with a straight line with its value and slope at the border.
        With `knot_step` <= 0, the average spectrum is smoothed with a spline
        of degree 5 with the default smoothing factor.

        `skip_windows` are wavelength intervals (Angstrom) with absorptions in
        the light of the calibration lamp: they are left out of the fit of the
        spline, and in the fiber flat of each fiber they are replaced by a
        straight line between both sides.

        A Savitzky-Golay filter (`window`, `degree`) of the average spectrum is
        only plotted in the intermediate results.
        """
        from scipy.signal import savgol_filter
        from scipy.interpolate import LSQUnivariateSpline, UnivariateSpline

        # Bad fibers
        fp_conf = FocalPlaneConf.from_img(rss_wl)
        good_idxs_mask = numpy.ones((fp_conf.nfibers,), dtype="bool")
        good_idxs_mask[[fibid - 1 for fibid in fp_conf.invalid_fibers()]] = False

        data0 = rss_wl[0].data
        # This extension was created by WLcalibrator
        wlmap = rss_wl["WLMAP"].data > 0

        # Level of each fiber
        level = numpy.median(data0[:, col1:col2], axis=1)
        lit = good_idxs_mask & (level > 0)
        ref_level = numpy.median(level[lit]) if lit.any() else 0.0
        valid_mask = lit & (level > min_level * ref_level)
        faint = numpy.flatnonzero(lit & ~valid_mask) + 1
        if faint.size:
            self.logger.warning("fibers too faint for the average spectrum: %s", faint.tolist())

        valid_level = level[valid_mask]
        covered = wlmap[valid_mask]
        with numpy.errstate(invalid="ignore", divide="ignore"):
            data_good = numpy.where(covered, data0[valid_mask] / valid_level[:, numpy.newaxis], 0.0)
        data_good[~numpy.isfinite(data_good)] = 0.0

        # number of fibers in the average of each wavelength
        mm = covered.sum(axis=0)
        mask_noinfo = mm < 1
        collapse = numpy.where(mask_noinfo, 1.0, data_good.sum(axis=0) / numpy.maximum(mm, 1))
        # Smoothing works bad very near the border (overshooting)
        collapse_smooth = savgol_filter(collapse, window, degree)
        collapse_smooth[mask_noinfo] = 1.0

        xx = numpy.arange(collapse.shape[0])
        skip = numpy.zeros(collapse.shape[0], dtype=bool)
        if skip_windows:
            hdr = rss_wl[0].header
            wl = hdr.get("CRVAL1", 0.0) + hdr.get("CDELT1", 1.0) * (xx + 1 - hdr.get("CRPIX1", 1.0))
            for w1, w2 in skip_windows:
                skip |= (wl >= w1) & (wl <= w2)
            self.logger.info("skipping %d pixels in %d absorption windows", skip.sum(), len(skip_windows))
        collapse_smooth_s = numpy.ones_like(collapse)
        if knot_step > 0:
            cdelt = abs(rss_wl[0].header.get("CDELT1", 1.0))
            coverage = mm / max(mm.max(), 1)
            fit = (coverage >= full_coverage) & ~skip
            x_fit = xx[fit]
            knots = _flat_knots(x_fit, knot_step / cdelt)
            self.logger.debug(
                "average spectrum smoothed with %d interior knots in pixels %d-%d", len(knots), x_fit[0], x_fit[-1]
            )
            spline = LSQUnivariateSpline(x_fit, collapse[fit], knots, k=3)
            collapse_smooth_s[~mask_noinfo] = _extend_linear(spline, xx[~mask_noinfo], x_fit[0], x_fit[-1])
            spline_label = f"spline, knots every {knot_step:g} A"
        else:
            spline = UnivariateSpline(xx[~mask_noinfo], collapse[~mask_noinfo], k=5)
            spline_label = "spline5"
            collapse_smooth_s[~mask_noinfo] = spline(xx[~mask_noinfo])

        if self.intermediate_results:
            from matplotlib.figure import Figure

            numpy.savetxt("collapse.txt", collapse)
            numpy.savetxt("mask_noinfo.txt", mask_noinfo)
            fig = Figure()
            ax = fig.subplots()
            ax.plot(xx, collapse, ".", label="collapsed")
            ax.plot(xx, collapse_smooth, "-", label=f"savgol{degree}")
            ax.plot(xx, collapse_smooth_s, "--", label=spline_label)
            ax.legend()
            fig.savefig("collapsed_smooth.png")

        # Divide each fiber in rss_wl by spectrum
        gmean = valid_level.mean()
        data1 = rss_wl[0].data / collapse_smooth_s
        data1 /= gmean
        if skip.any():
            interpolate_windows(data1, skip, wlmap)
        # Fill values with ones to avoid NaNs
        data2 = numpy.where(wlmap, data1, 1.0)

        self.logger.warning("Copy all extensions for the moment")
        rss_wl2 = fits.HDUList([hdu.copy() for hdu in rss_wl])
        rss_wl2[0].data = data2

        # quality flags of the fibers
        flags = flag_fibers(data2, wlmap, lit & ~valid_mask, good_idxs_mask, col1, col2, skip)
        fibers_hdr = rss_wl2["FIBERS"].header
        for fibid, fiber_flags in sorted(flags.items()):
            fibers_hdr[f"FIB{fibid:03d}QF"] = (",".join(fiber_flags), "Quality flags of the fiber flat")
            if not INVALID_FLAT_FLAGS.isdisjoint(fiber_flags):
                fibers_hdr[f"FIB{fibid:03d}_V"] = (False, "Fiber is invalid")
        deviant = [fibid for fibid, fiber_flags in flags.items() if "deviant" in fiber_flags]
        if deviant:
            self.logger.info("fibers with a flat different from their neighbors: %s", deviant)
        broken = [fibid for fibid, fiber_flags in flags.items() if "broken" in fiber_flags]
        if broken:
            self.logger.warning("fibers with a broken fiber flat, marked invalid: %s", broken)
        return rss_wl2

    def run_qc(self, recipe_input, recipe_result):
        """Quality of the fiber flat, from the flags of the fibers"""
        self.logger.info("start fiber flat recipe QC")
        with recipe_result.master_fiberflat.open() as hdul:
            recipe_result.qc = check_qc_fiberflat(hdul["FIBERS"].header)
        self.logger.info("Result QC is %s", recipe_result.qc)
        self.logger.info("end fiber flat recipe QC")
        return recipe_result

    def run(self, rinput):
        """Execute the recipe.

        Parameters
        ----------
        rinput : RecipeInput

        Returns
        -------
        RecipeResult

        """

        img = self.process_flat2d(rinput)
        self.save_intermediate_img(img, "reduced_image.fits")
        splitter1 = Splitter()
        calibrator_aper = ApertureExtractor(rinput.master_apertures, self.datamodel, offset=rinput.extraction_offset)
        splitter2 = Splitter()
        calibrator_wl = WavelengthCalibrator(rinput.master_wlcalib, self.datamodel)
        flipcor = FlipLR()

        img = splitter1(img)
        flat2d = splitter1.out  # Copy before extraction
        img = calibrator_aper(img)
        img = splitter2(img)
        rss_base = splitter2.out  # Copy before el calibration
        zero_order = []
        if rinput.remove_zero_order:
            traces = _trace_polynomials(rinput.master_apertures)
            zero_order = find_zero_order(flat2d[0].data, traces)
            for band in zero_order:
                cols = band["center"](numpy.array(band["rows"]))
                self.logger.info(
                    "zero order of the VPH at columns %.0f-%.0f, FWHM %.1f px, %.1f %% of the fibers between them",
                    cols.min(),
                    cols.max(),
                    band["fwhm"],
                    100 * band["amplitude"],
                )
            remove_zero_order(img[0].data, traces, zero_order)
        self.logger.debug("Flip RSS left-rigtht, before WL calibration")
        img = flipcor.run(img)
        # Calibrate in WL
        rss_wl = calibrator_wl(img)
        self.save_intermediate_img(rss_wl, "reduced_rss.fits")

        # Obtain flat field
        self.logger.info("Normalize flat field")
        skip_windows = None
        if rinput.skip_absorption_windows:
            resolution = str(rss_wl[0].header.get("VPH", ""))[:2]
            skip_windows = absorption_windows(resolution)
            if not skip_windows:
                self.logger.info("no absorption windows for VPH %s", rss_wl[0].header.get("VPH"))
        rss_wl2 = self.obtain_fiber_flat(
            rss_wl, window=rinput.smoothing_window, knot_step=rinput.smoothing_knot_step, skip_windows=skip_windows
        )
        rss_wl2[0].header = self.set_product_headers(rss_wl2[0].header)
        for band in zero_order:
            cols = band["center"](numpy.array(band["rows"]))
            rss_wl2[0].header["history"] = (
                f"Zero order of the VPH interpolated at columns {cols.min():.0f}-{cols.max():.0f} "
                f"(FWHM {band['fwhm']:.1f} px)"
            )
        result = self.create_result(master_fiberflat=rss_wl2, reduced_image=flat2d, reduced_rss=rss_base)
        return result

    def set_product_headers(self, hdr):
        """Set metadata in the FITS headers of the main product."""
        hdr = self.set_base_headers(hdr)
        hdr["NUMTYPE"] = ("MasterFiberFlat", "Product type")
        return hdr
