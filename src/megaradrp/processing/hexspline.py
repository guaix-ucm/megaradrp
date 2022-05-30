#
# Copyright 2017-2020 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""
Interpolation method based on:
'Hex-Splines: A Novel Spline Family for Hexagonal Lattices'
van de Ville et al. IEEE Transactions on Image Processing 2004, 13, 6
"""

from __future__ import division

import math
import numpy as np
from megaradrp.simulation.convolution import hex_c

# Hexagon constants
M_SQRT3 = math.sqrt(3)


def hexspline1(xx, yy):
    """Hexspline of order p=1"""
    return hex_c(xx, yy, rad=1.0 / M_SQRT3, ang=0.0)


def hexspline2(x1, x2):
    """Hexspline of order p=2"""

    # This function has been created with sympy

    x1 = np.asanyarray(x1)
    x2 = np.asanyarray(x2)
    Heaviside = np.heaviside

    M_SQRT3_D_3 = M_SQRT3 / 3
    M_SQRT3_D_2 = M_SQRT3 / 2
    M_SQRT3_D_6 = M_SQRT3 / 6
    M_SQRT3_2D3 = 2 * M_SQRT3_D_3
    ref = 1 / 2

    x_p0_m1 = x2 - 1
    x_p0_p1 = x2 + 1
    x_p0_m12 = x2 - 1 / 2
    x_p0_p12 = x2 + 1 / 2

    x_p0_scaled = M_SQRT3_D_3 * x2
    x_p0_m1_scaled = M_SQRT3_D_3 * x_p0_m1
    x_p0_p1_scaled = M_SQRT3_D_3 * x_p0_p1
    x_p0_m12_scaled = M_SQRT3_D_3 * x_p0_m12
    x_p0_p12_scaled = M_SQRT3_D_3 * x_p0_p12

    aux_p1 = x1 + x_p0_m1_scaled
    aux_p1_m3 = aux_p1 - M_SQRT3_D_3

    aux_m1 = x1 - x_p0_m1_scaled
    aux_m1_p3 = x1 - x_p0_m1_scaled + M_SQRT3_D_3

    aux_mm1 = x1 - x_p0_p1_scaled
    aux_mm1_m3 = aux_mm1 - M_SQRT3_D_3
    aux_pm1 = x1 + M_SQRT3 * (x_p0_p1) / 3
    aux_pm1_p3 = aux_pm1 + M_SQRT3_D_3

    aux_p12 = x1 + x_p0_m12_scaled
    aux_p12_m2 = aux_p12 - M_SQRT3_D_2
    aux_p12_p2 = aux_p12 + M_SQRT3_D_2
    aux_p12_p6 = aux_p12 + M_SQRT3_D_6
    aux_p12_m6 = aux_p12 - M_SQRT3_D_6

    aux_m12 = x1 - x_p0_m12_scaled
    aux_m12_m6 = aux_m12- M_SQRT3_D_6
    aux_m12_p6 = aux_m12+ M_SQRT3_D_6
    aux_m12_p2 = aux_m12+ M_SQRT3_D_2
    aux_m12_m2 = aux_m12- M_SQRT3_D_2

    aux_mm12 = x1 - x_p0_p12_scaled
    aux_mm12_p6 = aux_mm12 + M_SQRT3_D_6
    aux_mm12_m6 = aux_mm12 - M_SQRT3_D_6
    aux_mm12_m2 = aux_mm12 - M_SQRT3_D_2
    aux_mm12_p2 = aux_mm12 + M_SQRT3_D_2

    aux_pm12 = x1 + x_p0_p12_scaled
    aux_pm12_m6 = aux_pm12 - M_SQRT3_D_6
    aux_pm12_p6 = aux_pm12 + M_SQRT3_D_6
    aux_pm12_p2 = aux_pm12 + M_SQRT3_D_2
    aux_pm12_m2 = aux_pm12 - M_SQRT3_D_2

    aux_m0 = x1 - x_p0_scaled
    aux_m0_p23 = aux_m0 + M_SQRT3_2D3
    aux_m0_m23 = aux_m0 - M_SQRT3_2D3
    aux_m0_p3 = aux_m0 + M_SQRT3_D_3
    aux_m0_m3 = aux_m0 - M_SQRT3_D_3

    aux_p0 = x1 + x_p0_scaled
    aux_p0_m23 = aux_p0 - M_SQRT3_2D3
    aux_p0_p23 = aux_p0 + M_SQRT3_2D3
    aux_p0_p3 = aux_p0 + M_SQRT3_D_3
    aux_p0_m3 = aux_p0 - M_SQRT3_D_3

    h_p0 = Heaviside(x2, ref)
    h_p0_m1 = Heaviside(x_p0_m1, ref)
    h_p0_p1 = Heaviside(x_p0_p1, ref)
    h_p0_m12 = Heaviside(x_p0_m12, ref)
    h_p0_p12 = Heaviside(x_p0_p12, ref)

    h_m0 = Heaviside(-x2, ref)
    h_m0_p12 = Heaviside(1 / 2 - x2, ref)
    h_m0_m12 = Heaviside(-x2 - 1 / 2, ref)
    h_m0_p1 = Heaviside(1 - x2, ref)
    h_m0_m1 = Heaviside(-x2 - 1, ref)

    h_aux_m0 = Heaviside(aux_m0, ref)
    h_aux_m1 = Heaviside(aux_m1, ref)
    h_aux_p0 = Heaviside(aux_p0, ref)
    h_aux_p1 = Heaviside(aux_p1, ref)
    h_aux_m0_p23 = Heaviside(aux_m0_p23, ref)
    h_aux_m0_m23 = Heaviside(aux_m0_m23, ref)
    h_aux_p0_m23 = Heaviside(aux_p0_m23, ref)
    h_aux_p0_p23 = Heaviside(aux_p0_p23, ref)
    h_aux_m1_p3 = Heaviside(aux_m1_p3, ref)
    h_aux_p1_m3 = Heaviside(aux_p1_m3, ref)
    h_aux_m12_m6 = Heaviside(aux_m12_m6, ref)
    h_aux_m12_p2 = Heaviside(aux_m12_p2, ref)
    h_aux_p12_m2 = Heaviside(aux_p12_m2, ref)
    h_aux_p12_p6 = Heaviside(aux_p12_p6, ref)
    h_aux_mm12_m2 = Heaviside(aux_mm12_m2, ref)
    h_aux_mm12_p6 = Heaviside(aux_mm12_p6, ref)
    h_aux_pm12_m6 = Heaviside(aux_pm12_m6, ref)
    h_aux_pm12_p2 = Heaviside(aux_pm12_p2, ref)
    h_aux_mm1_m3 = Heaviside(aux_mm1_m3, ref)
    h_aux_pm1_p3 = Heaviside(aux_pm1_p3, ref)
    h_aux_mm1 = Heaviside(aux_mm1, ref)

    scale = M_SQRT3_2D3
    res = 4 * x2 * aux_m0 * h_p0 * h_aux_m0 - \
          4 * x2 * aux_p0 * h_m0 * h_aux_p0 + \
          x2 * aux_m0_m23 * h_p0 * h_aux_m0_m23 + \
          x2 * aux_m0_p23 * h_p0 * h_aux_m0_p23 - \
          x2 * aux_p0_m23 * h_m0 * h_aux_p0_m23 - \
          x2 * aux_p0_p23 * h_m0 * h_aux_p0_p23 + \
          x_p0_m1 * aux_m1_p3 * h_p0_m1 * h_aux_m1_p3 - \
          x_p0_m1 * aux_p1_m3 * h_m0_p1 * h_aux_p1_m3 - \
          2 * (x_p0_m12) * aux_m12_m6 * h_p0_m12 * h_aux_m12_m6 - \
          2 * (x_p0_m12) * aux_m12_p2 * h_p0_m12 * h_aux_m12_p2 + \
          2 * (x_p0_m12) * aux_p12_m2 * h_m0_p12 * h_aux_p12_m2 + \
          2 * (x_p0_m12) * aux_p12_p6 * h_m0_p12 * h_aux_p12_p6 - \
          2 * (x_p0_p12) * aux_mm12_m2 * h_p0_p12 * h_aux_mm12_m2 - \
          2 * (x_p0_p12) * aux_mm12_p6 * h_p0_p12 * h_aux_mm12_p6 + \
          2 * (x_p0_p12) * aux_pm12_m6 * h_m0_m12 * h_aux_pm12_m6 + \
          2 * (x_p0_p12) * aux_pm12_p2 * h_m0_m12 * h_aux_pm12_p2 + \
          (x_p0_p1) * (aux_mm1_m3) * h_p0_p1 * h_aux_mm1_m3 - \
          (x_p0_p1) * (aux_pm1_p3) * h_m0_m1 * h_aux_pm1_p3 + \
          M_SQRT3 * ((aux_m0) ** 2 * h_p0 * h_aux_m0 +
                     (aux_p0) ** 2 * h_m0 * h_aux_p0 +
                     (aux_m1) ** 2 * h_aux_m1 * h_p0_m1 / 2 +
                     (aux_p1) ** 2 * h_m0_p1 * h_aux_p1 / 2 +
                     (aux_mm1) ** 2 * h_aux_mm1 * h_p0_p1 / 2 +
                     (aux_mm1) ** 2 * h_aux_mm1 * h_m0_m1 / 2 +
                     (aux_m0_m23) ** 2 * h_p0 * h_aux_m0_m23 / 2 +
                     (aux_m0_m3) ** 2 * h_p0 * Heaviside(aux_m0_m3, ref) / 2 +
                     (aux_m0_p3) ** 2 * h_p0 * Heaviside(aux_m0_p3, ref) / 2 +
                     (aux_m0_p23) ** 2 * h_p0 * h_aux_m0_p23 / 2 +
                     (aux_p0_m23) ** 2 * h_m0 * h_aux_p0_m23 / 2 +
                     (aux_p0_m3) ** 2 * h_m0 * Heaviside(aux_p0_m3, ref) / 2 +
                     (aux_p0_p3) ** 2 * h_m0 * Heaviside(aux_p0_p3, ref) / 2 +
                     (aux_p0_p23) ** 2 * h_m0 * h_aux_p0_p23 / 2 -
                     (aux_m12_m2) ** 2 * h_p0_m12 * Heaviside(aux_m12_m2, ref) / 2 -
                     (aux_m12_m6) ** 2 * h_p0_m12 * h_aux_m12_m6 / 2 -
                     (aux_m12_p6) ** 2 * h_p0_m12 * Heaviside(aux_m12_p6, ref) / 2 -
                     (aux_m12_p2) ** 2 * h_p0_m12 * h_aux_m12_p2 / 2 -
                     (aux_p12_m2) ** 2 * h_m0_p12 * h_aux_p12_m2 / 2 -
                     (aux_p12_m6) ** 2 * h_m0_p12 * Heaviside(aux_p12_m6, ref) / 2 -
                     (aux_p12_p6) ** 2 * h_m0_p12 * h_aux_p12_p6 / 2 -
                     (aux_p12_p2) ** 2 * h_m0_p12 * Heaviside(aux_p12_p2, ref) / 2 -
                     (aux_mm12_m2) ** 2 * h_p0_p12 * h_aux_mm12_m2 / 2 -
                     (aux_mm12_m6) ** 2 * h_p0_p12 * Heaviside(aux_mm12_m6, ref) / 2 -
                     (aux_mm12_p6) ** 2 * h_p0_p12 * h_aux_mm12_p6 / 2 -
                     (aux_mm12_p2) ** 2 * h_p0_p12 * Heaviside(aux_mm12_p2, ref) / 2 -
                     (aux_pm12_m2) ** 2 * h_m0_m12 * Heaviside(aux_pm12_m2, ref) / 2 -
                     (aux_pm12_m6) ** 2 * h_m0_m12 * h_aux_pm12_m6 / 2 -
                     (aux_pm12_p6) ** 2 * h_m0_m12 * Heaviside(aux_pm12_p6, ref) / 2 -
                     (aux_pm12_p2) ** 2 * h_m0_m12 * h_aux_pm12_p2 / 2)
    res = res * scale
    return res


def _hexspline2_t1(x1, x2):
    """Hexspline of order p=2 in region t1"""

    # This function has been created with sympy

    x1 = np.asanyarray(np.abs(x1))
    x2 = np.asanyarray(np.abs(x2))

    M_SQRT3_D_3 = M_SQRT3 / 3
    M_SQRT3_D_2 = M_SQRT3 / 2
    M_SQRT3_D_6 = M_SQRT3 / 6
    M_SQRT3_2D3 = 2 * M_SQRT3_D_3

    x_p0_m12 = x2 - 1 / 2
    x_p0_p12 = x2 + 1 / 2

    x_p0_scaled = M_SQRT3_D_3 * x2

    x_p0_m12_scaled = M_SQRT3_D_3 * x_p0_m12
    x_p0_p12_scaled = M_SQRT3_D_3 * x_p0_p12

    aux_p12 = x1 + x_p0_m12_scaled
    aux_p12_p2 = aux_p12 + M_SQRT3_D_2
    aux_p12_p6 = aux_p12 + M_SQRT3_D_6
    aux_mm12 = x1 - x_p0_p12_scaled
    aux_mm12_p2 = aux_mm12 + M_SQRT3_D_2
    aux_m0 = x1 - x_p0_scaled
    aux_m0_p23 = aux_m0 + M_SQRT3_2D3
    aux_m0_p3 = aux_m0 + M_SQRT3_D_3

    scale = M_SQRT3_2D3
    res = x2 * aux_m0_p23  + \
          2 * x_p0_m12 * aux_p12_p6  + \
          M_SQRT3 * (
                     (aux_m0_p3) ** 2 / 2 +
                     (aux_m0_p23) ** 2 / 2 -
                     (aux_p12_p6) ** 2 / 2 -
                     (aux_p12_p2) ** 2  / 2 -
                     (aux_mm12_p2) ** 2 / 2
                  )
    res = res * scale
    return res


def _hexspline2_regions(x1, x2):
    """Hexspline regions of order p=2"""
    x1 = np.asanyarray(x1)
    x2 = np.asanyarray(x2)
    Heaviside = np.heaviside

    M_SQRT3_D_3 = M_SQRT3 / 3
    M_SQRT3_D_2 = M_SQRT3 / 2
    M_SQRT3_D_6 = M_SQRT3 / 6
    M_SQRT3_2D3 = 2 * M_SQRT3_D_3
    ref = 1 / 2

    x_p0_m1 = x2 - 1
    x_p0_p1 = x2 + 1
    x_p0_m12 = x2 - 1 / 2
    x_p0_p12 = x2 + 1 / 2

    x_p0_scaled = M_SQRT3_D_3 * x2
    x_p0_m1_scaled = M_SQRT3_D_3 * x_p0_m1
    x_p0_p1_scaled = M_SQRT3_D_3 * x_p0_p1
    x_p0_m12_scaled = M_SQRT3_D_3 * x_p0_m12
    x_p0_p12_scaled = M_SQRT3_D_3 * x_p0_p12

    aux_p1 = x1 + x_p0_m1_scaled
    aux_p1_m3 = aux_p1 - M_SQRT3_D_3

    aux_m1 = x1 - x_p0_m1_scaled
    aux_m1_p3 = x1 - x_p0_m1_scaled + M_SQRT3_D_3

    aux_mm1 = x1 - x_p0_p1_scaled
    aux_mm1_m3 = aux_mm1 - M_SQRT3_D_3
    aux_pm1 = x1 + M_SQRT3 * (x_p0_p1) / 3
    aux_pm1_p3 = aux_pm1 + M_SQRT3_D_3

    aux_p12 = x1 + x_p0_m12_scaled
    aux_p12_m2 = aux_p12 - M_SQRT3_D_2
    aux_p12_p2 = aux_p12 + M_SQRT3_D_2
    aux_p12_p6 = aux_p12 + M_SQRT3_D_6
    aux_p12_m6 = aux_p12 - M_SQRT3_D_6

    aux_m12 = x1 - x_p0_m12_scaled
    aux_m12_m6 = aux_m12- M_SQRT3_D_6
    aux_m12_p6 = aux_m12+ M_SQRT3_D_6
    aux_m12_p2 = aux_m12+ M_SQRT3_D_2
    aux_m12_m2 = aux_m12- M_SQRT3_D_2

    aux_mm12 = x1 - x_p0_p12_scaled
    aux_mm12_p6 = aux_mm12 + M_SQRT3_D_6
    aux_mm12_m6 = aux_mm12 - M_SQRT3_D_6
    aux_mm12_m2 = aux_mm12 - M_SQRT3_D_2
    aux_mm12_p2 = aux_mm12 + M_SQRT3_D_2

    aux_pm12 = x1 + x_p0_p12_scaled
    aux_pm12_m6 = aux_pm12 - M_SQRT3_D_6
    aux_pm12_p6 = aux_pm12 + M_SQRT3_D_6
    aux_pm12_p2 = aux_pm12 + M_SQRT3_D_2
    aux_pm12_m2 = aux_pm12 - M_SQRT3_D_2

    aux_m0 = x1 - x_p0_scaled
    aux_m0_p23 = aux_m0 + M_SQRT3_2D3
    aux_m0_m23 = aux_m0 - M_SQRT3_2D3
    aux_m0_p3 = aux_m0 + M_SQRT3_D_3
    aux_m0_m3 = aux_m0 - M_SQRT3_D_3

    aux_p0 = x1 + x_p0_scaled
    aux_p0_m23 = aux_p0 - M_SQRT3_2D3
    aux_p0_p23 = aux_p0 + M_SQRT3_2D3
    aux_p0_p3 = aux_p0 + M_SQRT3_D_3
    aux_p0_m3 = aux_p0 - M_SQRT3_D_3

    r = {}
    h_p0 = Heaviside(x2, ref)
    h_p0_m1 = Heaviside(x_p0_m1, ref)
    h_p0_p1 = Heaviside(x_p0_p1, ref)
    h_p0_m12 = Heaviside(x_p0_m12, ref)
    h_p0_p12 = Heaviside(x_p0_p12, ref)

    h_m0 = Heaviside(-x2, ref)
    h_m0_p12 = Heaviside(1 / 2 - x2, ref)
    h_m0_m12 = Heaviside(-x2 - 1 / 2, ref)
    h_m0_p1 = Heaviside(1 - x2, ref)
    h_m0_m1 = Heaviside(-x2 - 1, ref)

    h_aux_m0 = Heaviside(aux_m0, ref)
    h_aux_m1 = Heaviside(aux_m1, ref)
    h_aux_p0 = Heaviside(aux_p0, ref)
    h_aux_p1 = Heaviside(aux_p1, ref)
    h_aux_m0_p23 = Heaviside(aux_m0_p23, ref)
    h_aux_m0_m23 = Heaviside(aux_m0_m23, ref)
    h_aux_p0_m23 = Heaviside(aux_p0_m23, ref)
    h_aux_p0_p23 = Heaviside(aux_p0_p23, ref)
    h_aux_m1_p3 = Heaviside(aux_m1_p3, ref)
    h_aux_p1_m3 = Heaviside(aux_p1_m3, ref)
    h_aux_m12_m6 = Heaviside(aux_m12_m6, ref)
    h_aux_m12_p2 = Heaviside(aux_m12_p2, ref)
    h_aux_p12_m2 = Heaviside(aux_p12_m2, ref)
    h_aux_p12_p6 = Heaviside(aux_p12_p6, ref)
    h_aux_mm12_m2 = Heaviside(aux_mm12_m2, ref)
    h_aux_mm12_p6 = Heaviside(aux_mm12_p6, ref)
    h_aux_pm12_m6 = Heaviside(aux_pm12_m6, ref)
    h_aux_pm12_p2 = Heaviside(aux_pm12_p2, ref)
    h_aux_mm1_m3 = Heaviside(aux_mm1_m3, ref)
    h_aux_pm1_p3 = Heaviside(aux_pm1_p3, ref)
    h_aux_mm1 = Heaviside(aux_mm1, ref)

    h_aux_m0_m3 = Heaviside(aux_m0_m3, ref)
    h_aux_m0_p3 = Heaviside(aux_m0_p3, ref)
    h_aux_p0_m3 = Heaviside(aux_p0_m3, ref)
    h_aux_p0_p3 = Heaviside(aux_p0_p3, ref)

    h_aux_m12_m2 = Heaviside(aux_m12_m2, ref)
    h_aux_m12_p6 = Heaviside(aux_m12_p6, ref)
    h_aux_p12_m6 = Heaviside(aux_p12_m6, ref)
    h_aux_p12_p2 = Heaviside(aux_p12_p2, ref)
    h_aux_mm12_m6 = Heaviside(aux_mm12_m6, ref)
    h_aux_mm12_p2 = Heaviside(aux_mm12_p2, ref)
    h_aux_pm12_m2 = Heaviside(aux_pm12_m2, ref)
    h_aux_pm12_p6 = Heaviside(aux_pm12_p6, ref)

    inter = locals().copy()
    for key in inter:
        if key.startswith('h_'):
            r[key] = inter[key]
    return r


def hexspline_support(xx, yy, p):
    """Support of hexspline of order p

    The support is an hexagon of area Omega p^2
    Omega = M_SQRT3 / 2
    """
    return hex_c(xx, yy, rad=p / M_SQRT3, ang=0.0)


def hexspline_bbox(p):
    """Bounding box of hexspline of order p

    The support is an hexagon of area Omega p^2
    Omega = M_SQRT3 / 2
    """
    # (x1, x2), (y1, y2)
    # TODO: review numbers
    x1 = p / M_SQRT3
    y1 = p / 2.0
    return (-x1, x1), (-y1, y1)


def hexspline_gauss(xx, yy, p):
    """Approximation of the p-order hexspline by a bivariate normal

    https://miplab.epfl.ch/pub/vandeville0202.pdf

    Least-squares spline resampling to a hexagonal lattice
    Signal Processing:Image Communication17 (2002) 393-408

    Eq B.4
    """
    from scipy.stats import multivariate_normal
    zz = np.array([xx, yy])
    zz = np.moveaxis(zz, 0, -1)
    in_support = hexspline_support(xx, yy, p)

    coeff = 5 * M_SQRT3 / 144
    omega = M_SQRT3 / 2
    sigma = coeff * np.eye(2)
    covar = p * sigma / omega

    out = omega * multivariate_normal.pdf(zz, mean=[0, 0], cov=covar)
    out[~in_support] = 0
    return out


def OsincH_2pi(x, y):
    t1 = np.cos(-x / (2 * M_SQRT3) + y / 2) - np.cos(x / M_SQRT3)
    t2 = np.cos(x / (2 * M_SQRT3) + y / 2) - np.cos(x / M_SQRT3)
    ts = t1 / (x + M_SQRT3 * y) + t2 / (x - M_SQRT3 * y)
    return 2 * M_SQRT3 * ts / x


def sincH(x, y):
    """Generalization of the sinc function to a hexagonal grid

    The function is 1 in (0,0) and 0 in all the knots of the grid
    """
    # This expression seems numerically unstable near 0
    # it should be 1
    x = np.asanyarray(x)
    y = np.asanyarray(y)

    xx = 2 * np.pi * np.where(x == 0, 1.0e-20, x)
    yy = 2 * np.pi * np.where(y == 0, 1.0e-20, y)

    t3 = np.cos(xx / M_SQRT3)
    t1 = np.cos(-xx / (2 * M_SQRT3) + yy / 2.0) - t3
    t2 = np.cos(xx / (2 * M_SQRT3) + yy / 2.0) - t3
    ts = t1 / (xx + M_SQRT3 * yy) + t2 / (xx - M_SQRT3 * yy)
    return 4 * ts / xx


# Convolution, compute kernel
def rescaling_kernel(p, scale=1):
    """Rescaling kernel from hexgrid to rectangular grid"""
    from megaradrp.simulation.convolution import setup_grid
    from scipy import signal
    from scipy.interpolate import RectBivariateSpline

    Dx = 0.005
    Dy = 0.005
    DA = Dx * Dy
    # TODO: the support should be computed from the scale and p
    xsize = ysize = 3.0
    xx, yy, xs, ys, xl, yl = setup_grid(xsize, ysize, Dx, Dy)

    detR0 = M_SQRT3 / 2
    detR1 = scale * scale

    # index of bspline
    n = p - 1

    rect_kernel = signal.bspline(xx / scale, n) * signal.bspline(yy / scale, n) / detR1

    hex1 = hexspline1(xx, yy)

    if p == 1:
        hex_kernel = hex1
    elif p == 2:
        hex2 = signal.fftconvolve(hex1, hex1, mode='same') * DA / detR0
        hex_kernel = hex2
    elif p == 3:
        hex2 = signal.fftconvolve(hex1, hex1, mode='same') * DA / detR0
        hex3 = signal.fftconvolve(hex2, hex1, mode='same') * DA / detR0
        hex_kernel = hex3
    else:
        raise ValueError('p>3 not implemented')

    kernel = signal.fftconvolve(rect_kernel, hex_kernel, mode='same') * DA
    rbs = RectBivariateSpline(xs, ys, kernel)
    return rbs
