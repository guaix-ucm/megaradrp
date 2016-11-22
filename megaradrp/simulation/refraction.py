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

    PS = -10474.0 + 116.43 * T - 0.43284 * T ** 2 + 0.00053840 * T ** 3

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
