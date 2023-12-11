#
# Copyright 2015-2023 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

import sys
if sys.version_info[:2] <= (3, 8):
    from typing import Sequence  # for typing
else:
    from collections.abc import Sequence  # for typing

from numina.modeling.gaussbox import GaussBox
import numpy as np

from .base import ModelDescription


class GaussBoxModelDescription(ModelDescription):

    def __init__(self, fixed_center=True, init_simple=False, npix=5, sigma=3.0):
        super().__init__("gaussbox", GaussBox, fixed_center=fixed_center,
                         params_fit=['amplitude', 'mean', 'stddev'],
                         params_save=['mean', 'stddev'],
                         deg_save=[3, 5])
        self.npix = npix
        self.sigma = sigma
        self.init_simple = init_simple

    @classmethod
    def create_from_centers(cls, fixed_center=True, **kwargs):

        # Parameters of __init__ for GaussBox
        init_simple = kwargs.get('init_simple', False)
        npix = kwargs.get('npix', 5)
        sigma = kwargs.get('sigma', 3.0)

        # TODO: this method can be completely general
        # using inspect to obtain the signature of __init__

        # obj = cls.__new__(cls)
        # obj.__init__(fixed_center=fixed_center,
        #     init_simple=init_simple,
        #     npix=npix, sigma=sigma)

        obj = GaussBoxModelDescription(
            fixed_center=fixed_center, init_simple=init_simple,
            npix=npix, sigma=sigma
        )

        return obj

    def _init_1(self, column: np.array, centers: Sequence[float]):
        """init simple"""
        npix = self.npix
        ecenters = np.ceil(centers - 0.5).astype('int')
        nfib = len(ecenters)

        offset = npix // 2
        yfit = np.empty((npix, nfib))
        xfit = np.arange(npix) - offset

        for i in range(npix):
            yfit[i, :] = column[ecenters + (i - offset)]

        coeff = np.polyfit(xfit, np.log(yfit), deg=2)
        c, b, a = coeff

        sig2 = -1 / (2 * c)
        mu = b / sig2
        ampl = np.exp(a - mu ** 2)

        params = {
            'mean': centers,
            'stddev': np.sqrt(sig2),
            'amplitude': ampl,
            'hpix': 0.5 * np.ones_like(centers),
        }
        return params

    def _init_2(self, column: np.array, centers: Sequence[float]) -> dict:
        ecenters = np.ceil(centers - 0.5).astype('int')
        nfib = len(centers)
        sigma = self.sigma
        params = {
            'mean': centers,
            'stddev': sigma * np.ones_like(centers),
            'amplitude': [column[ecenters[i]] / 0.25 for i in range(nfib)],
            'hpix': 0.5 * np.ones_like(centers),
        }
        return params

    def init_values(self, column: np.array, centers: Sequence[float]):
        if self.init_simple:
            return self._init_1(column, centers)
        else:
            return self._init_2(column, centers)

    def params_fixed(self, values, fibid, col):
        result = {
            'hpix': True
        }
        if self.fixed_center:
            result['mean'] = True
        return result

    def params_bounds(self, values, fibid, col):

        isigma = values['stddev']
        imean = values['mean']

        p_bounds = {
            "stddev": (isigma - 0.5, isigma + 0.5),
            "amplitude": (0, None)
        }
        if not self.fixed_center:
            p_bounds['mean'] = (imean - 0.5, imean + 0.5)
        return p_bounds

    def fiber_center(self, values: dict) -> float:
        center = values['mean']
        return center
