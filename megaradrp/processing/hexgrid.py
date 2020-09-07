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

from __future__ import print_function

import math

import numpy as np

# Normalized hexagon geometry
M_SQRT3 = math.sqrt(3)
H_HEX = 0.5
R_HEX = 1 / M_SQRT3
A_HEX = 0.5 * H_HEX * R_HEX
HA_HEX = 6 * A_HEX  # detR0 == ha


def calc_matrix(nrow, ncol, grid_type=2):
    """

    Parameters
    ----------
    nrow : int
    ncol : int
    grid_type : int

    Returns
    -------

    """

    rr0 = np.array([[M_SQRT3 / 2, 0], [-0.5, 1]])  # Unit scale

    if grid_type == 2:
        f = 0
    else:
        f = 1

    kcol = []
    krow = []
    for i in range(ncol):
        s = (i + f * (i % 2)) // 2
        for j in range(nrow):
            kcol.append(i)
            krow.append(j+s)

    sl = np.array([kcol, krow])  # x y
    r0l = np.dot(rr0, sl)
    # r0l = R0 @ sl
    return r0l


def calc_lcb_grid(scale=1.0):
    """

    Parameters
    ----------
    scale : float

    Returns
    -------

    """

    grid_type = 2  # Values for MEGARA
    ncol = 27  #
    nrow = 21  #
    r0l = calc_matrix(nrow, ncol, grid_type=grid_type)
    # r0l = R0 @ sl
    spos_x = scale * (r0l[0] - r0l[0].max() / 2)
    spos_y = scale * (r0l[1] - r0l[1].max() / 2)

    return spos_x, spos_y


def hexgrid_extremes(r0l, target_scale):
    """

    Parameters
    ----------
    r0l
    target_scale : float

    Returns
    -------

    """
    # geometry
    # ha_hex = 6 * a_hex # detR0 == ha
    # compute extremes of hexgrid to rectangular grid
    # with pixel size 'scale'
    x0min, y0min = r0l.min(axis=1)
    x0max, y0max = r0l.max(axis=1)
    y1min = y0min - H_HEX
    y1max = y0max + H_HEX
    x1min = x0min - R_HEX
    x1max = x0max + R_HEX

    j1min = int(math.floor(x1min / target_scale + 0.5))
    i1min = int(math.floor(y1min / target_scale + 0.5))
    j1max = int(math.ceil(x1max / target_scale - 0.5))
    i1max = int(math.ceil(y1max / target_scale - 0.5))
    return (i1min, i1max), (j1min, j1max)


def connected6(i, j):
    """Connected pixels in a hexagonal grid"""
    for a, b in [(0, 1), (0, -1), (1, 0), (-1, 0), (1, -1), (-1, 1)]:
        u, v = i + a, j + b
        yield u, v


def neigh6(i, j, level=1):
    """neighborhood in 6-connection"""
    # square
    for m in range(0, level + 1):
        for n in range(0, level + 1):
            yield i - m, j + n
    for k in range(level):
        for ll in range(1, level + 1):
            yield i + 1 + k, j - 1 - k + ll
            yield i + 1 + k - ll, j - 1 - k
        yield i + 1 + k, j - 1 - k
