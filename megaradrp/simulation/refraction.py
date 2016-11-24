
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

"""Computing differential atmospheric refraction"""

from __future__ import print_function

import astropy.units as u

import math

import numpy

#D2R = math.pi / 180.0

@u.quantity_input(zenith_distance=u.deg)
def differential(
        zenith_distance,
        wl,
        wl_reference,
        temperature,
        relative_humidity,
        pressure
):
    # Z is zenith distance
    # lambda0 reference wavelength [um]

    #TC = 11.5                        #Temperature [C]
    # RH = 14.5                        #Relative Humidity [%]
    # P = 743.0                        #Pressure [mbar]

    ZD = zenith_distance.to(u.rad).value
    T = temperature.to(u.K, equivalencies=u.temperature()).value
    Lambda = wl.to(u.micron).value
    Lambda0 = wl_reference.to(u.micron).value
    RH = relative_humidity / 100
    P = pressure.to(u.bar).value * 1e3 # milibar

    return _dar_dimensionless_(ZD, T, RH, P, Lambda, Lambda0) * u.rad


def _dar_dimensionless_(ZD, T, RH, P, wl, lambda_ref):

    # print(ZD, T, RH, P, wl, Lambda0)

    PS = -10474.0 + T * (116.43 - T *(0.43284 + 0.00053840 * T))

    P2 = RH * PS
    P1 = P - P2
    # print(P1)
    D1 = P1 / T * (1.0 + P1 * (57.90 * 1.0E-8 - (9.3250 * 1.0E-4 / T) + (0.25844 / T ** 2)))
    D2 = P2 / T * (1.0 + P2 * (1.0 + 3.7E-4 * P2) * (-2.37321E-3 + (2.23366 / T) - (710.792 / T ** 2) + (7.75141E4 / T ** 3)))

    N0_1 = _poly1(D1, D2, 1.0 / lambda_ref)

    N_1 = _poly1(D1, D2, 1.0 / wl)

    DR = math.tan(ZD) * (N0_1 - N_1)
    return DR


def _poly1(D1, D2, S):

    N_1 = 1.0E-8 * ((2371.34 + 683939.7 / (130 - S ** 2) + 4547.3 / (38.9 - S ** 2)) * D1 + (
        6487.31 + 58.058 * S ** 2 - 0.71150 * S ** 4 + 0.08851 * S ** 6) * D2)
    return N_1


if __name__ == '__main__':

    wl = numpy.linspace(0.3, 1.0, num=15) * u.micron

    #wl = 0.6 * u.micron

    print(wl)
    print(differential(60*u.deg, wl, 0.5 * u.micron, 11.5 * u.deg_C, 44.2, 772.2e-3 * u.bar).to(u.arcsec))

# [ 0.3   0.35  0.4   0.45  0.5   0.55  0.6   0.65  0.7   0.75  0.8   0.85
#   0.9   0.95  1.  ] micron
# [ 139.21843488   75.42767344   39.3957048    16.18259071    0.
#   -11.86954952  -20.8895519   -27.92780533  -33.53521547  -38.07938831
#  -41.81507414  -44.92415491  -47.53967765  -49.76096379  -51.66345799] arcsec