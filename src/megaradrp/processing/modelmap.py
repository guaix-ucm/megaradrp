#
# Copyright 2015-2025 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

import bisect
import collections
import logging
import math
import warnings

from collections.abc import Sequence  # for typing

from astropy.modeling import fitting
from astropy.modeling.functional_models import Const1D
from numina.util.objimport import import_object
import numpy

from .modeldesc.base import ModelDescription
from megaradrp.processing.modeldesc import config

_logger = logging.getLogger(__name__)


def calc1d_model(
    model_desc: ModelDescription,
    column: numpy.ndarray,
    centers: Sequence[float],
    valid: Sequence[int],
    col: int,
    lateral=2,
    reject=3,
    nloop=1,
    fit_background=False,
    width_limits=None,
    smooth_widths=0,
) -> dict:
    """Fit a sum of profiles along a 1D column vector

    Each fiber is fitted together with its `lateral` neighbours at each side,
    after removing the profiles of the fibers up to `reject` positions away,
    in `nloop` loops over the fibers in random order.

    Parameters
    ----------
    fit_background : bool
        Fit a constant background in each local fit, instead of a background
        fixed at zero. Where the profiles overlap, as in defocused zones of
        the detector, the light between the fibers otherwise goes into the
        widths and amplitudes of the profiles.
    width_limits : tuple of float, optional
        Absolute limits of the width (stddev) of the profiles. Without them,
        the width is bounded around its current value in each loop, so the
        bounds move from loop to loop.
    smooth_widths : int
        If positive, after the loops the width of each fiber is replaced by
        a smooth value from the fibers around it on the slit (see
        smooth_profile_widths), and a last loop fits the amplitudes with
        those widths fixed. Neighbouring fibers that overlap otherwise trade
        light, alternating wide and narrow profiles.
    """

    nfib = len(centers)
    if nfib != len(valid):
        raise ValueError("len(valid) fibers must equal to len(centers)")

    model_cls = model_desc.model_cls

    # shape of LCB image
    # nfibers = 623
    nrows = 4112

    xl = numpy.arange(nrows, dtype=float)

    if column.shape != xl.shape:
        raise ValueError(f"boxd1d must have len {nrows}")

    # 'current' contains the values of the parameters
    # for each fiber
    # it is initially filled with initial values
    current = model_desc.init_values_per_profile(column, centers, valid)

    # total fits is 2 * lateral + 1
    total = reject + lateral
    npix = 6
    loops = nloop + (1 if smooth_widths > 0 else 0)
    for il in range(loops):
        fixed_width = smooth_widths > 0 and il == loops - 1
        if fixed_width:
            smooth_profile_widths(current, valid, smooth_widths)
        # print('loop', il, datetime.datetime.now())
        # permutation of the valid fibers
        permutated_fibers = numpy.random.permutation(valid)
        for fib_id_reorder in permutated_fibers:
            #
            # print(f'val is {val}')
            # Compute the center of the fiber from the fitted parameters
            _logger.debug("fib %s", fib_id_reorder)
            params = current[fib_id_reorder]
            _logger.debug("params of this fiber %s", params)
            center = model_desc.fiber_center(params)
            _logger.debug("center of this fiber %s", center)
            # Region containing the fibers and its neighbours
            m0 = int(math.ceil(center - 0.5))
            m1 = max(0, m0 - npix * total)
            m2 = m0 + npix * total

            _logger.debug("m1:m2 is %s:%s", m1, m2)
            # Cut regions in the arrays
            yt = column[m1:m2].copy()
            xt = xl[m1:m2]

            pos = bisect.bisect(valid, fib_id_reorder)
            _logger.debug("pos bisect %s", pos)
            candidates = valid[max(pos - total - 1, 0) : pos + total]
            _logger.debug("candidates around my fiber %s", list(candidates))
            dis_f = [abs(c - fib_id_reorder) for c in candidates]
            _logger.debug("distances around my fibers %s", dis_f)
            # for p in range(max(0, val-S), min(nfib, val+S+1)):
            fit_ids = []
            for cfib, fib_df in zip(candidates, dis_f):
                if fib_df > lateral:
                    # remove contribution
                    yt -= model_cls.evaluate(xt, **current[cfib])
                else:
                    fit_ids.append(cfib)
            _logger.debug("candidates around my fiber %s", fit_ids)
            # Extract values to create initials
            # We add profiles to a constant value of zero

            model = Const1D(amplitude=0.0)
            model.amplitude.fixed = not fit_background

            # For each fiber to be fitted
            for fib_id in fit_ids:
                # The values of the parameters
                permutated_fibers = current[fib_id]
                # Parameters to be fixed
                p_fixed = model_desc.params_fixed(permutated_fibers, fib_id, col)
                # Bounded parameters
                p_bounds = model_desc.params_bounds(permutated_fibers, fib_id, col)
                if width_limits is not None:
                    p_bounds = dict(p_bounds, stddev=tuple(width_limits))
                if fixed_width:
                    p_fixed = dict(p_fixed, stddev=True)
                # The values of the parameters in a tuple
                args = tuple(permutated_fibers[name] for name in model_cls.param_names)
                # The model for this fiber
                newg = model_cls(*args, fixed=p_fixed, bounds=p_bounds)
                # Add to the current model
                model = model + newg

            # Fitter for this loop
            fitter = fitting.LevMarLSQFitter()
            # Fit the model, the model is new in each loop, it is not copied
            model_fitted = fitter(model, xt, yt, inplace=True)

            # Extract the parameters of the model and
            # update 'current' with them
            for val_idx, fibid in enumerate(fit_ids, 1):
                # This will overwrite fits with fits performed later
                # if fibid == val:
                fvalue = {}
                for name in model_cls.param_names:
                    fvalue[name] = getattr(model_fitted, f"{name}_{val_idx}").value

                for name in model_cls.param_names:
                    current[fibid][name] = fvalue[name]

    return current


def smooth_profile_widths(current, valid, half):
    """Replace the width of each fiber by a smooth value from the fibers around it on the slit

    The width of each fiber is first averaged with that of the next one,
    which removes an alternation of wide and narrow profiles, and then the
    median of 2 * half + 1 of those averages around each fiber is taken.
    `current` maps fiber ids to the parameters of their profiles, changed in
    place.
    """
    order = sorted(valid)
    widths = numpy.array([current[fibid]["stddev"] for fibid in order], dtype=float)
    if widths.size < 2:
        return
    pairs = 0.5 * (widths[:-1] + widths[1:])
    # one value per fiber: the average with the next fiber, the last one with the previous
    pairs = numpy.append(pairs, pairs[-1])
    for k, fibid in enumerate(order):
        current[fibid]["stddev"] = float(numpy.median(pairs[max(k - half, 0) : k + half + 1]))


def calc_matrix(wshape, model, params, valid, clip=1.0e-6, extra=10):
    from scipy.sparse import csr_matrix

    # calc w
    # this is valid for 1 column (there are more broadcasts for N)
    g_mean = params["mean"]

    begpix = numpy.ceil(g_mean - 0.5).astype("int")

    steps = numpy.arange(-extra, extra)
    # ref is a bidimensional matrix +-10 pixels around the trace
    ref = begpix + steps[:, numpy.newaxis]

    for key in model.param_names:
        if key not in params:
            model_param = getattr(model, key)
            params[key] = model_param.value

    rr = model.evaluate(ref, **params)
    rrb = begpix - extra

    # This was sending warnings. Is there a NaN somewhere?
    with numpy.errstate(invalid="ignore"):
        rr[rr < clip] = 0.0

    # Calc Ws matrix, the profile of each valid fiber in its column,
    # starting in its row rrb, without the zero values
    block = rr.shape[0]
    idx = numpy.asarray(valid) - 1
    rows = rrb[idx] + numpy.arange(block)[:, numpy.newaxis]
    cols = numpy.broadcast_to(idx, rows.shape)
    values = rr[:, idx]
    nonzero = values != 0
    wcol = csr_matrix((values[nonzero], (rows[nonzero], cols[nonzero])), shape=wshape)
    return wcol


def calc_matrix_adapt(wshape, col, model, params, valid, clip=1.0e-6, extra=10):
    """Adapted function to return the column number"""

    # For parallel processing
    wcol = calc_matrix(wshape, model, params, valid, clip=clip, extra=extra)

    return col, wcol


def calc_matrix_cols(model_map, datashape, processes=0):

    dnrow, dncol = datashape
    nfibs = model_map.total_fibers

    shape = (nfibs, dncol)
    wshape = (dnrow, nfibs)
    #
    params_reorder = collections.defaultdict(lambda: numpy.zeros(shape))
    model_name = "undefined"
    xcol = numpy.arange(dncol)
    for fibermodel in model_map.contents:
        if fibermodel.valid:
            row = fibermodel.fibid - 1
            interpol_params = fibermodel.model["params"]
            # it only makes sense to use one model specification
            model_name = fibermodel.model["model_name"]
            for name, value in interpol_params.items():
                params_reorder[name][row] = value(xcol)

    if model_name not in config:
        raise ValueError(f"model name {model_name} is not defined")

    objpath = config[model_name]
    model_class = import_object(objpath)
    modeldesc = model_class()
    model = modeldesc.model_cls

    # shift center parameter according to global_offset
    mask = numpy.zeros((model_map.total_fibers,), dtype="bool")
    valid_r = [(f.fibid - 1) for f in model_map.contents if f.valid]
    mask[valid_r] = True

    params_reorder_c = modeldesc.fiber_center(params_reorder)
    mean_at_ref = params_reorder_c[mask, model_map.ref_column]
    offset = model_map.global_offset(mean_at_ref)
    params_reorder_c[mask, :] = params_reorder_c[mask, :] + offset[:, numpy.newaxis]

    wcols = {}
    valid = [f.fibid for f in model_map.contents if f.valid]

    # parameters per colum
    params_col = {}
    for col in xcol:
        mpar = {}
        for key in model.param_names:
            if key not in params_reorder:
                model_param = getattr(model, key)
                mpar[key] = model_param.value
            else:
                mpar[key] = params_reorder[key][:, col]
        params_col[col] = mpar

    if processes < 2:
        for col in xcol:
            result = calc_matrix(wshape, model, params_col[col], valid, clip=1e-6, extra=10)
            wcols[col] = result
    else:
        from numina.util.parallel import logging_pool

        # the results are obtained inside the block, the pool
        # is terminated at the end
        with logging_pool(processes) as pool:
            results = [
                pool.apply_async(
                    calc_matrix_adapt,
                    args=(wshape, col, model, params_col[col], valid),
                    kwds={"clip": 1e-6, "extra": 10},
                )
                for col in xcol
            ]

            for p in results:
                col, wcol = p.get()
                wcols[col] = wcol

    return wcols


def aper_extract(model_map, wcols, img):
    from scipy.sparse.linalg import lsqr

    n0 = model_map.total_fibers
    n1 = img.shape[1]
    rss = numpy.zeros((n0, n1))
    for key, val in wcols.items():
        yl = img[:, key]
        res = lsqr(val, yl)
        rss[:, key] = res[0]

    return rss


# flags of the quality control of the model map
MODEL_FLAG_CROSSING = "crossing"
MODEL_FLAG_NOT_PARALLEL = "not_parallel"
MODEL_FLAG_WIDTH = "width_not_smooth"
MODEL_FLAG_AMPLITUDE = "amplitude_not_smooth"


def _robust_limit(values, nsigma, floor):
    """median + nsigma robust sigmas of the finite `values`, and at least `floor`"""
    finite = values[numpy.isfinite(values)]
    if finite.size == 0:
        return floor
    med = numpy.median(finite)
    sigma = 1.4826 * numpy.median(numpy.abs(finite - med))
    return max(med + nsigma * sigma, floor)


def _neighbour_prediction(values, boxes, use=None):
    """Prediction of each row of `values` from its neighbours of the same block

    The neighbours are the closest fibers of the same block at each side
    that are in `use` (all by default). The prediction is their mean, or the
    only neighbour for the fibers at the border of a block (an extrapolation
    from two neighbours doubles the noise); NaN for a fiber without
    neighbours.
    """
    nfib = values.shape[0]
    use = numpy.ones(nfib, dtype=bool) if use is None else use
    pred = numpy.full(values.shape, numpy.nan)
    for i in range(nfib):
        prev = next((j for j in range(i - 1, -1, -1) if boxes[j] == boxes[i] and use[j]), None)
        succ = next((j for j in range(i + 1, nfib) if boxes[j] == boxes[i] and use[j]), None)
        if prev is not None and succ is not None:
            pred[i] = 0.5 * (values[prev] + values[succ])
        elif prev is not None:
            pred[i] = values[prev]
        elif succ is not None:
            pred[i] = values[succ]
    return pred


def _flag_outliers(values, boxes, metric, nsigma, floor, niter=3):
    """Fibers whose `metric` against the prediction from their neighbours is an outlier

    The flagged fibers are not used as neighbours in the next iteration, so
    that a bad fiber does not flag its neighbours.
    """
    use = numpy.ones(values.shape[0], dtype=bool)
    bad = numpy.zeros(values.shape[0], dtype=bool)
    for _ in range(niter):
        measure = metric(values, _neighbour_prediction(values, boxes, use))
        limit = _robust_limit(measure[use], nsigma, floor)
        new_bad = measure > limit
        if numpy.array_equal(new_bad, bad):
            break
        bad = new_bad
        use = ~bad
    return numpy.flatnonzero(bad)


def _parallel_metric(centres, pred):
    """Range along the detector of the centre minus the prediction"""
    return numpy.ptp(centres - pred, axis=1)


def _width_metric(widths, pred):
    """Largest difference of the width from the prediction"""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return numpy.nanmax(numpy.abs(widths - pred), axis=1)


def _amplitude_metric(flux, pred):
    """Largest local relative change of the ratio between the fiber and the prediction

    The ratio removes the spectrum of the lamp; its median over 21 columns is
    compared with its median over 401 columns.
    """
    from scipy.ndimage import median_filter

    with numpy.errstate(invalid="ignore", divide="ignore"):
        ratio = flux / pred
    rough = numpy.full(flux.shape[0], numpy.nan)
    for i, r in enumerate(ratio):
        if not numpy.all(numpy.isfinite(r)) or numpy.median(r) <= 0:
            continue
        local = median_filter(r, 21, mode="nearest")
        broad = median_filter(r, 401, mode="nearest")
        rough[i] = numpy.max(numpy.abs(local / broad - 1))
    return rough


def check_model_map(
    model_map,
    rss=None,
    cols=None,
    nsigma=5.0,
    min_gap=0.5,
    parallel_floor=0.5,
    width_floor=0.3,
    amplitude_floor=0.05,
    bad_fraction=0.05,
    logger=None,
):
    """Quality control of the geometry of a model map

    The checks, for each valid fiber, compare its model with those of its
    neighbours of the same block of the pseudo-slit:

    - the traces do not cross: the gap between neighbouring centres is
      larger than `min_gap` times the median gap, at all the columns;
    - each trace is approximately parallel to its neighbours: the range,
      along the detector, of its centre minus the prediction from its
      neighbours;
    - the width varies smoothly: the largest difference from the width
      predicted by its neighbours;
    - the amplitude varies smoothly along the fiber: with the RSS of the flat
      extracted with the model (`rss`), the ratio between the fiber and the
      prediction from its neighbours (which removes the spectrum of the lamp)
      is compared with a smooth version of itself. The amplitude is not
      compared between fibers: a fiber can be fainter than its neighbours.

    A fiber fails a check if its value is an outlier (beyond the median plus
    `nsigma` robust sigmas of all the fibers) and larger than the floor of the
    check (in pixels, or relative for the amplitude). It gets a flag in
    `flags`, but its model is not changed. The quality of the model map is
    GOOD without flagged fibers, PARTIAL with some, and BAD if they are more
    than `bad_fraction` of the valid fibers.

    Returns
    -------
    dict
        The fibids of the fibers with each flag
    """
    from numina.types.qc import QC

    logger = logger or _logger
    models = [m for m in model_map.contents if m.valid]
    if not models:
        model_map.quality_control = QC.BAD
        return {}
    if cols is None:
        cols = numpy.arange(200, 3901, 100)
    cols = numpy.asarray(cols, dtype=float)
    centres = numpy.array([m.model["params"]["mean"](cols) for m in models])
    order = numpy.argsort(centres[:, len(cols) // 2])
    models = [models[i] for i in order]
    centres = centres[order]
    fibids = numpy.array([m.fibid for m in models])
    boxes = numpy.array([m.boxid for m in models])
    found = {}

    def flag(name, indices):
        for i in indices:
            models[i].add_flag(name)
            found.setdefault(name, set()).add(int(fibids[i]))

    # crossings
    gaps = numpy.diff(centres, axis=0)
    same_box = boxes[1:] == boxes[:-1]
    typical = numpy.median(gaps[same_box]) if numpy.any(same_box) else numpy.median(gaps)
    close = numpy.flatnonzero(gaps.min(axis=1) < min_gap * typical)
    flag(MODEL_FLAG_CROSSING, set(close) | set(close + 1))

    # parallelism
    flag(MODEL_FLAG_NOT_PARALLEL, _flag_outliers(centres, boxes, _parallel_metric, nsigma, parallel_floor))

    # width
    if all("stddev" in m.model["params"] for m in models):
        widths = numpy.array([m.model["params"]["stddev"](cols) for m in models])
        flag(MODEL_FLAG_WIDTH, _flag_outliers(widths, boxes, _width_metric, nsigma, width_floor))

    # amplitude along each fiber, relative to its neighbours
    if rss is not None:
        flux = numpy.asarray(rss, dtype=float)[fibids - 1][:, int(cols[0]) : int(cols[-1])]
        flag(MODEL_FLAG_AMPLITUDE, _flag_outliers(flux, boxes, _amplitude_metric, nsigma, amplitude_floor))

    flagged = set().union(*found.values()) if found else set()
    for name, values in sorted(found.items()):
        logger.warning("model map: fibers with flag %r: %d %s", name, len(values), sorted(values))
    if len(flagged) > bad_fraction * len(models):
        model_map.quality_control = QC.BAD
    elif flagged:
        model_map.quality_control = QC.PARTIAL
    else:
        model_map.quality_control = QC.GOOD
    return {name: sorted(values) for name, values in found.items()}
