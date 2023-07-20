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


from astropy.modeling.functional_models import Moffat1D
from numina.modeling.gaussbox import GaussBox
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
            raise ValueError('deg_save and params_save must have the same length')

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
            'amplitude': ampl
        }
        return params

    def _init_2(self, column: np.array, centers: Sequence[float]) -> dict:
        ecenters = np.ceil(centers - 0.5).astype('int')
        nfib = len(centers)
        sigma = self.sigma
        params = {
            'mean': centers,
            'stddev': sigma * np.ones_like(centers),
            'amplitude': [column[ecenters[i]] / 0.25 for i in range(nfib)]
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
