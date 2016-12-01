#
# Copyright 2016 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# Megara DRP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Megara DRP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Megara DRP.  If not, see <http://www.gnu.org/licenses/>.
#

import math
from astropy.modeling.functional_models import Gaussian2D, Moffat2D


class AtmosphereModel(object):

    def __init__(self, twilight, nightsky, seeing, extinction, refraction):
        self.tw_interp = twilight
        self.ng_interp = nightsky
        self.ext_interp = extinction
        self.seeing = seeing
        self.refraction_model = refraction

    def twilight_spectrum(self, wl_in):
        """Twilight spectrum"""
        return self.tw_interp(wl_in)

    def night_spectrum(self, wl_in):
        """Night spectrum"""
        return self.ng_interp(wl_in)

    def extinction(self, wl_in):
        """Night extinction"""
        return self.ext_interp(wl_in)

    def refraction(self, z, wl, ref):
        return self.refraction_model.refraction(z, wl, ref)


class SeeingSizeModel(object):
    def __init__(self, wl, r0):
        self._r0 = r0
        self._wl0 = wl

    def fwhm(self, wl, zd):
        return 1.2 * wl / self.r0(wl, zd)

    def r0(self, wl, zd):
        return self._r0 * (wl / self._wl0) ** 1.2 * math.cos(zd)

    def profile(self, fwhm):
        return generate_gaussian_profile(fwhm)


class ConstSeeing(object):
    def __init__(self, seeing):
        self._s = seeing

    def fwhm(self, wl, zd):
        return self._s

    def profile(self, fwhm):
        return generate_gaussian_profile(fwhm)


def generate_gaussian_profile(seeing_fwhm):
    """Generate a normalized Gaussian profile from its FWHM"""
    FWHM_G = 2 * math.sqrt(2 * math.log(2))
    sigma = seeing_fwhm / FWHM_G
    amplitude = 1.0 / (2 * math.pi * sigma * sigma)
    seeing_model = Gaussian2D(amplitude=amplitude,
                              x_mean=0.0,
                              y_mean=0.0,
                              x_stddev=sigma,
                              y_stddev=sigma)
    return seeing_model


def generate_moffat_profile(seeing_fwhm, alpha):
    """Generate a normalized Moffat profile from its FWHM and alpha"""

    scale = 2 * math.sqrt(2**(1.0 / alpha) - 1)
    gamma = seeing_fwhm / scale
    amplitude = 1.0 / math.pi * (alpha - 1) / gamma**2
    seeing_model = Moffat2D(amplitude=amplitude,
                            x_mean=0.0,
                            y_mean=0.0,
                            gamma=gamma,
                            alpha=alpha)
    return seeing_model


def generate_lorentz_profile(seeing_fwhm):
    """Generate a normalized Lorent profile from its FWHM"""

    return generate_moffat_profile(seeing_fwhm, alpha=1.5)
