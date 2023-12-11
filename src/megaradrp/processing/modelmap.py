#
# Copyright 2015-2023 Universidad Complutense de Madrid
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
import sys

if sys.version_info[:2] <= (3, 8):
    from typing import Sequence  # for typing
else:
    from collections.abc import Sequence  # for typing

from astropy.modeling import fitting
from astropy.modeling.functional_models import Const1D
from numina.util.objimport import import_object
import numpy

from .modeldesc.base import ModelDescription
from megaradrp.processing.modeldesc import config


def calc1d_model(model_desc: ModelDescription, column: numpy.ndarray,
                 centers: Sequence[float],
                 valid: Sequence[int], col: int,
                 lateral=2, reject=3, nloop=1) -> dict:
    """Fit a sum of profiles along a 1D column vector"""

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
    for il in range(nloop):
        # print('loop', il, datetime.datetime.now())
        # permutation of the valid fibers
        permutated_fibers = numpy.random.permutation(valid)
        for fib_id_reorder in permutated_fibers:
            #
            # print(f'val is {val}')
            # Compute the center of the fiber from the fitted parameters
            logging.debug(f'fib {fib_id_reorder}')
            params = current[fib_id_reorder]
            logging.debug(f"params of this fiber {params}")
            center = model_desc.fiber_center(params)
            logging.debug(f"center of this fiber {center}")
            # Region containing the fibers and its neighbours
            m0 = int(math.ceil(center - 0.5))
            m1 = max(0, m0 - npix * total)
            m2 = m0 + npix * total

            logging.debug(f'm1:m2 is {m1}:{m2}')
            # Cut regions in the arrays
            yt = column[m1:m2].copy()
            xt = xl[m1:m2]

            pos = bisect.bisect(valid, fib_id_reorder)
            logging.debug(f'pos bisect {pos}')
            candidates = valid[max(pos - total - 1, 0): pos + total]
            logging.debug(f'candidates around my fiber  {list(candidates)}')
            dis_f = [abs(c - fib_id_reorder) for c in candidates]
            logging.debug(f'distances around my fibers {dis_f}')
            # for p in range(max(0, val-S), min(nfib, val+S+1)):
            fit_ids = []
            for cfib, fib_df in zip(candidates, dis_f):
                if fib_df > lateral:
                    # remove contribution
                    yt -= model_cls.evaluate(xt, **current[cfib])
                else:
                    fit_ids.append(cfib)
            logging.debug(f'candidates around my fiber {fit_ids}')
            # Extract values to create initials
            # We add profiles to a constant value of zero

            model = Const1D(amplitude=0.0)
            model.amplitude.fixed = True

            # For each fiber to be fitted
            for fib_id in fit_ids:
                # The values of the parameters
                permutated_fibers = current[fib_id]
                # Parameters to be fixed
                p_fixed = model_desc.params_fixed(
                    permutated_fibers, fib_id, col)
                # Bounded parameters
                p_bounds = model_desc.params_bounds(
                    permutated_fibers, fib_id, col)
                # The values of the parameters in a tuple
                args = tuple(permutated_fibers[name]
                             for name in model_cls.param_names)
                # The model for this fiber
                newg = model_cls(*args, fixed=p_fixed, bounds=p_bounds)
                # Add to the current model
                model = model + newg

            # Fitter for this loop
            fitter = fitting.LevMarLSQFitter()
            # Fit the model
            model_fitted = fitter(model, xt, yt)

            # Extract the parameters of the model and
            # update 'current' with them
            for val_idx, fibid in enumerate(fit_ids, 1):
                # This will overwrite fits with fits performed later
                # if fibid == val:
                fvalue = {}
                for name in model_cls.param_names:
                    fvalue[name] = getattr(
                        model_fitted, f'{name}_{val_idx}').value

                for name in model_cls.param_names:
                    current[fibid][name] = fvalue[name]

    return current


def calc_matrix(wshape, model, params, valid, clip=1.0e-6, extra=10):
    from scipy.sparse import lil_matrix

    # calc w
    # this is valid for 1 column (there are more broadcasts for N)
    g_mean = params['mean']

    begpix = numpy.ceil(g_mean - 0.5).astype('int')

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
    with numpy.errstate(invalid='ignore'):
        rr[rr < clip] = 0.0

    # Calc Ws matrix
    block, valid_nfib = rr.shape

    w_init = lil_matrix(wshape)
    for fibid in valid:
        idx = fibid - 1
        w_init[rrb[idx]:rrb[idx] + block, idx] = rr[:, idx, numpy.newaxis]

    wcol = w_init.tocsr()
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
    model_name = 'undefined'
    xcol = numpy.arange(dncol)
    for fibermodel in model_map.contents:
        if fibermodel.valid:
            row = fibermodel.fibid - 1
            interpol_params = fibermodel.model['params']
            # it only makes sense to use one model specification
            model_name = fibermodel.model['model_name']
            for name, value in interpol_params.items():
                params_reorder[name][row] = value(xcol)

    if model_name not in config:
        raise ValueError(f"model name {model_name} is not defined")

    objpath = config[model_name]
    model_class = import_object(objpath)
    modeldesc = model_class()
    model = modeldesc.model_cls

    # shift center parameter according to global_offset
    mask = numpy.zeros((model_map.total_fibers,), dtype='bool')
    valid_r = [(f.fibid - 1) for f in model_map.contents if f.valid]
    mask[valid_r] = True

    params_reorder_c = modeldesc.fiber_center(params_reorder)
    mean_at_ref = params_reorder_c[mask, model_map.ref_column]
    offset = model_map.global_offset(mean_at_ref)
    params_reorder_c[mask, :] = params_reorder_c[mask, :] + \
        offset[:, numpy.newaxis]

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
            result = calc_matrix(wshape, model, params_col[col], valid, clip=1e-6,
                                 extra=10)
            wcols[col] = result
    else:
        import multiprocessing as mp
        pool = mp.Pool(processes=processes)

        results = [pool.apply_async(
            calc_matrix_adapt,
            args=(wshape, col, model, params_col[col], valid),
            kwds={'clip': 1e-6, 'extra': 10}
        ) for col in xcol]

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
