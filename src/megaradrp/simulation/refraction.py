#
# Copyright 2016-2017 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Computing differential atmospheric refraction"""

from __future__ import print_function

import astropy.units as u

import math

import numpy


@u.quantity_input(zenith_distance=u.micron)
def ref_index0(wl):
    wl0 = wl.to(u.micron).value
    wl1 = (1.0 / wl0)**2
    return 64.328 + 29498.1 / (146 - wl1) + 255.4 / (41 - wl1)


@u.quantity_input(zenith_distance=u.micron)
def ref_index1(wl, p=760, t=15):
    # 1 mmHg = 133.322387415
    n6_0 = ref_index0(wl)
    factor1 = p * (1 + (1.049 - 0.0157 * t) * 1e-6 * p)
    factor2 = 720.883 * (1 + 0.003661 * t)
    n6_1 = n6_0 * factor1 / factor2
    return n6_1


@u.quantity_input(zenith_distance=u.micron)
def ref_index2(wl, p=760, t=15, f=0.0):
    n6_1 = ref_index1(wl, p, t)
    wl0 = wl.to(u.micron).value
    wl1 = (1.0 / wl0) ** 2
    cont1 = 0.0624 - 0.000680 * wl1
    cont2 = 720.883 * (1 + 0.003661 * t)
    n6_2 = n6_1 - cont1 / cont2 * f
    return n6_2


@u.quantity_input(zenith_distance=u.deg)
def differential_p(
        zenith_distance,
        wl,
        wl_reference,
        temperature,
        pressure,
        relative_humidity,
):
    """Differential refraction as given by 1982PASP...94..715F """

    zd = zenith_distance.to(u.rad).value
    p = pressure.to(u.Pascal).value / 133.322387415
    t = temperature.to(u.deg_C, equivalencies=u.temperature()).value
    f = relative_humidity * p
    nl = 1 + 1e-6 * ref_index2(wl, p, t, f)
    n0 = 1 + 1e-6 * ref_index2(wl_reference, p, t, f)

    delt_r = (nl - n0) * math.tan(zd)

    return delt_r * u.rad


class DifferentialRefractionModel(object):
    def __init__(self, temperature, pressure, relative_humidity):
        self.temperature = temperature
        self.pressure = pressure
        self.relative_humidity = relative_humidity

    def refraction(self, z, wl, ref):
        return differential_p(z, wl, ref, self.temperature, self.pressure, self.relative_humidity)


if __name__ == '__main__':
    from astropy.units import cds
    cds.enable()

    wl = numpy.linspace(0.3, 1.0, num=15) * u.micron

    #wl = 0.6 * u.micron

    print(wl)
    zdis = 60*u.deg
    rel = 8.0 / 600.0
    ref = 0.5 * u.micron
    temp = 11.5 * u.deg_C
    press = 600 * cds.mmHg
    print(differential_p(zdis, wl, ref, temp, press, rel).to(u.arcsec))


# [ 0.3   0.35  0.4   0.45  0.5   0.55  0.6   0.65  0.7   0.75  0.8   0.85
#   0.9   0.95  1.  ] micron
# [ 139.21843488   75.42767344   39.3957048    16.18259071    0.
#   -11.86954952  -20.8895519   -27.92780533  -33.53521547  -38.07938831
#  -41.81507414  -44.92415491  -47.53967765  -49.76096379  -51.66345799] arcsec