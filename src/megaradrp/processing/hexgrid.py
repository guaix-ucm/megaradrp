#
# Copyright 2017-2025 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

"""
Interpolation method based on:
'Hex-Splines: A Novel Spline Family for Hexagonal Lattices'
van de Ville et al. IEEE Transactions on Image Processing 2004, 13, 6
"""

import math

import numpy as np
from scipy.spatial import cKDTree

# Normalized hexagon geometry
M_SQRT3 = math.sqrt(3)
H_HEX = 0.5
R_HEX = 1 / M_SQRT3
A_HEX = 0.5 * H_HEX * R_HEX
HA_HEX = 6 * A_HEX  # detR0 == ha
# Offsets (x, y) to the 6 neighbors of a point in the grid of unit spacing
NEIGH_OFFSETS = np.array(
    [
        [0.0, M_SQRT3 / 2, M_SQRT3 / 2, 0.0, -M_SQRT3 / 2, -M_SQRT3 / 2],
        [1.0, 0.5, -0.5, -1.0, -0.5, 0.5],
    ]
)


def calc_matrix(nrow, ncol, grid_type=2):
    """
    Positions of the points of a hexagonal grid of unit spacing

    Parameters
    ----------
    nrow : int
        Number of points in each column
    ncol : int
        Number of columns
    grid_type : int, optional
        Type of grid, 2 is the grid of MEGARA, where the odd columns
        are shifted down by half a point

    Returns
    -------
    numpy.ndarray
        Positions (x, y) of the points, an array of shape (2, nrow * ncol)

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
            krow.append(j + s)

    sl = np.array([kcol, krow])  # x y
    r0l = np.dot(rr0, sl)
    # r0l = R0 @ sl
    return r0l


def calc_lcb_grid(scale=1.0):
    """
    Positions of the spaxels of the LCB, centered in the origin

    Parameters
    ----------
    scale : float, optional
        Distance between spaxels

    Returns
    -------
    tuple of numpy.ndarray
        Positions x and y of the 567 spaxels

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
    Extremes of the rectangular grid that covers a hexagonal grid

    The rectangular grid has square pixels of size `target_scale`, the
    center of the pixel (i, j) is (j * target_scale, i * target_scale).
    The pixels cover the hexagons around all the points.

    Parameters
    ----------
    r0l : numpy.ndarray
        Positions (x, y) of the points of the hexagonal grid, an array
        of shape (2, npoints), in units of the distance between points
    target_scale : float
        Size of the pixels of the rectangular grid, in units of the
        distance between points

    Returns
    -------
    tuple
        ((i1min, i1max), (j1min, j1max)), the minimum and maximum indices
        of the rows (y) and of the columns (x) of the rectangular grid

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


def hexgrid_outer_points(r0l):
    """
    Points of the hexagonal grid next to a set of points, but not in the set

    These are the points around the border of the set, and the
    points of the holes inside it.

    Parameters
    ----------
    r0l : numpy.ndarray
        Positions (x, y) of the points of the set, an array of shape
        (2, npoints), in units of the distance between points

    Returns
    -------
    outer : numpy.ndarray
        Positions (x, y) of the outer points, an array of shape
        (2, nouter)
    neighbors : list of list of int
        For each outer point, the indices in `r0l` of its neighbors

    """
    tree = cKDTree(r0l.T)
    cand = (r0l[:, :, np.newaxis] + NEIGH_OFFSETS[:, np.newaxis, :]).reshape(2, -1)
    dist, _ = tree.query(cand.T)
    cand = cand[:, dist > 0.5]
    # a point is a candidate once for each of its neighbors
    groups = cKDTree(cand.T).query_ball_point(cand.T, 0.5)
    keep = [idx for idx, group in enumerate(groups) if idx == min(group)]
    outer = cand[:, keep]
    neighbors = [list(nb) for nb in tree.query_ball_point(outer.T, 1.5)]
    return outer, neighbors


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
