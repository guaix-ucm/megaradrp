#
# Copyright 2015-2023 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import abc
from collections.abc import Sequence  # for typing
import bisect
import math

from astropy.modeling import fitting
from astropy.modeling.functional_models import Const1D, Moffat1D
from numina.modeling.gaussbox import GaussBox, gauss_box_model
import numpy as np
#import numpy.typing as npt


class ModelDescription(metaclass=abc.ABCMeta):

    def __init__(self, name, model_cls, fixed_center=True):
        self.model_cls = model_cls
        self.nfib = 1
        self.name = name
        self.fixed_center = fixed_center

    @abc.abstractmethod
    def init_values(self, column) -> dict:
        """Compute initial values for the parameters of the model"""
        raise NotImplementedError

    def init_values_per_profile(self, column, fibers: Sequence[int]) -> dict:
        """Rearrange the parameters of the profiles.

        Convert a dictionary with Np key parameters, with Nf values each
        into a dictionary with Nf keys and Np parameters each.
        """
        vals1 = self.init_values(column)

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


class MoffatModelDescription(ModelDescription):

    def __init__(self, centers1d, fixed_center=True):
        super().__init__("moffat", Moffat1D, fixed_center=fixed_center)
        self.nfib = len(centers1d)
        self.centers1d = centers1d

    def init_values(self, column):
        _params = {
            'x_0': self.centers1d,
            'gamma': 3.9 * np.ones_like(self.centers1d),
            'amplitude': 1 * np.ones_like(self.centers1d),
            'alpha': 2.9 * np.ones_like(self.centers1d)
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


class GaussBoxModelDescription(ModelDescription):

    def __init__(self, centers1d, fixed_center=True, init_simple=False, npix=5, sigma=3.0):
        super().__init__("gaussbox", GaussBox, fixed_center=fixed_center)
        self.nfib = len(centers1d)
        self.centers1d = centers1d
        self.npix = npix
        self.sigma = sigma
        self.init_simple = init_simple

    def _init_1(self, column):
        """init simple"""
        npix = self.npix
        ecenters = np.ceil(self.centers1d - 0.5).astype('int')
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
            'mean': self.centers1d,
            'stddev': np.sqrt(sig2),
            'amplitude': ampl
        }
        return params

    def _init_2(self, column):
        ecenters = np.ceil(self.centers1d - 0.5).astype('int')
        nfib = len(self.centers1d)
        sigma = self.sigma
        params = {
            'mean': self.centers1d,
            'stddev': sigma * np.ones_like(self.centers1d),
            'amplitude': [column[ecenters[i]] / 0.25 for i in range(nfib)]
        }
        return params

    def init_value(self, column):
        if self.init_simple:
            return self._init_1(column)
        else:
            return self._init_2(column)

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
