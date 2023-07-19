#
# Copyright 2015-2023 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import bisect
from collections.abc import Sequence  # for typing
import logging
import math

from astropy.modeling import fitting
from astropy.modeling.functional_models import Const1D
import numpy as np

from .modeldesc import ModelDescription


def calc1d_M(model_desc: ModelDescription, boxd1d: np.ndarray, valid: Sequence[int], col: int,
             lateral=2, reject=3, nloop=1) -> dict:
    """Fit a sum of profiles along a 1D column vector"""

    nfib = model_desc.nfib
    if nfib != len(valid):
        raise ValueError("len(valid) fibers must equal to len(centers)")

    model_cls = model_desc.model_cls

    # shape of LCB image
    nfibers = 623
    nrows = 4112

    xl = np.arange(nrows, dtype=float)

    if boxd1d.shape != xl.shape:
        raise ValueError(f"boxd1d must have len {nrows}")

    # 'current' contains the values of the parameters
    # for each fiber
    # it is initially filled with initial values
    current = model_desc.init_values_per_profile(boxd1d, valid)

    # total fits is 2 * lateral + 1
    total = reject + lateral
    npix = 6
    for il in range(nloop):
        # print('loop', il, datetime.datetime.now())
        # permutation of the valid fibers
        permutated_fibers = np.random.permutation(valid)
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
            yt = boxd1d[m1:m2].copy()
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
                p_fixed = model_desc.params_fixed(permutated_fibers, fib_id, col)
                # Bounded parameters
                p_bounds = model_desc.params_bounds(permutated_fibers, fib_id, col)
                # The values of the parameters in a tuple
                args = tuple(permutated_fibers[name] for name in model_cls.param_names)
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
                    fvalue[name] = getattr(model_fitted, f'{name}_{val_idx}').value

                for name in model_cls.param_names:
                    current[fibid][name] = fvalue[name]

    return current
