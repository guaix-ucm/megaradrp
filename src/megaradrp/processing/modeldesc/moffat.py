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

from astropy.modeling.functional_models import Moffat1D
import numpy as np

from .base import ModelDescription


class MoffatModelDescription(ModelDescription):

    def __init__(self, fixed_center=True):
        super().__init__("moffat", Moffat1D, fixed_center=fixed_center)

    def init_values(self, column: np.array, centers: Sequence[float]) -> dict:
        _params = {
            'x_0': centers,
            'gamma': 3.9 * np.ones_like(centers),
            'amplitude': 1 * np.ones_like(centers),
            'alpha': 2.9 * np.ones_like(centers)
        }
        return _params

    def params_fixed(self, values, fibid, col):
        result = {
            'alpha': False
        }
        if self.fixed_center:
            result['x_0'] = True
        return result

    def params_bounds(self, values, fibid, col):

        igamma = values['gamma']
        imean = values['x_0']

        p_bounds = {
            "gamma": (igamma - 0.5, igamma + 0.5),
            "amplitude": (0, None)
        }
        if not self.fixed_center:
            p_bounds['x_0'] = (imean - 0.5, imean + 0.5)
        return p_bounds

    def fiber_center(self, values: dict) -> float:
        center = values['x_0']
        return center
