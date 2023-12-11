#
# Copyright 2015-2023 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

import abc
import sys

if sys.version_info[:2] <= (3, 8):
    from typing import Sequence  # for typing
else:
    from collections.abc import Sequence  # for typing

import numpy as np


class ModelDescription(metaclass=abc.ABCMeta):

    def __init__(self, name, model_cls, fixed_center=True, params_fit=None, params_save=None, deg_save=None):
        self.model_cls = model_cls
        self.nfib = 1
        self.name = name
        self.fixed_center = fixed_center

        if params_fit is None:
            self.params_fit = self.model_cls.param_names
        else:
            self._check_params(params_fit)
            self.params_fit = params_fit

        if params_save is None:
            self.params_save = self.params_fit
        else:
            self._check_params(params_save)
            self.params_save = params_save
        if deg_save is None:
            # default, cubic spline for everything
            self.deg_save = [3 for _ in self.params_save]
        else:
            self.deg_save = deg_save

        if len(self.deg_save) != len(self.params_save):
            raise ValueError(
                'deg_save and params_save must have the same length')

    def _check_params(self, param_collection):
        """Check that names in param_c are valid param names"""
        for p in param_collection:
            if p not in self.model_cls.param_names:
                raise ValueError(f"parameter {p} not in model parameters")

    @abc.abstractmethod
    def init_values(self, column: np.array, centers: Sequence[float]) -> dict:
        """Compute initial values for the parameters of the model"""
        raise NotImplementedError

    def init_values_per_profile(self, column: np.array, centers: Sequence[float], fibers: Sequence[int]) -> dict:
        """Rearrange the parameters of the profiles.

        Convert a dictionary with Np key parameters, with Nf values each
        into a dictionary with Nf keys and Np parameters each.
        """
        vals1 = self.init_values(column, centers)

        current = {}
        for idx, fibid in enumerate(fibers):
            cf = {}
            for param_name in self.model_cls.param_names:
                cf[param_name] = vals1[param_name][idx]
            current[fibid] = cf

        return current

    @abc.abstractmethod
    def params_fixed(self, values, fibid, col):
        raise NotImplementedError

    @abc.abstractmethod
    def params_bounds(self, values, fibid, col):
        raise NotImplementedError

    @abc.abstractmethod
    def fiber_center(self, values: dict) -> float:
        raise NotImplementedError
