#
# Copyright 2015-2018 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import numpy as np
import math

from astropy.modeling.functional_models import Fittable2DModel, Parameter


_HEX_SCALE = 0.25 * math.sqrt(3.0)


class HexagonA(Fittable2DModel):

    amplitude = Parameter(default=1)
    x_0 = Parameter(default=0)
    y_0 = Parameter(default=0)
    radius = Parameter(default=1)
    angle = Parameter(default=0)

    @staticmethod
    def evaluate(x, y, amplitude, x_0, y_0, radius, angle):
        # Recenter
        xo = x - x_0
        yo = y - y_0
        # Rotate
        cs = math.cos(angle)
        ss = math.sin(angle)
        xr = cs * xo + ss * yo
        yr = -ss * xo + cs * yo
        # Scale
        d = radius * 2
        dx = np.abs(xr) / d
        dy = np.abs(yr) / d
        a = 0.25 * math.sqrt(3.0)
        sel1 = (dy <= a)
        sel2 = (a * dx + 0.25 * dy <= 0.5 * a)
        return np.select([np.logical_and(sel1, sel2)],
                         [amplitude], 0)

    @property
    def bounding_box(self):
        a = 0.25 * math.sqrt(3.0)
        d = 2 * self.radius
        return ((self.y_0 - a * d, self.y_0 + a * d),
                (self.x_0 - self.radius, self.x_0 + self.radius))


# Returns True in inside hexagon, rectangle and square
def hex_c(x, y, rad, ang=0.0):
    # Rotate coordinates
    cs = math.cos(ang)
    ss = math.sin(ang)
    xr = cs * x + ss *y
    yr = -ss * x + cs *y
    d = rad * 2
    dx = np.abs(xr) / d
    dy = np.abs(yr) / d
    a = 0.25 * math.sqrt(3.0)
    return (dy <= a) & (a*dx + 0.25*dy <= 0.5*a)


def rect_c(x, y, lx, ly):
    dx = np.abs(x) / lx
    dy = np.abs(y) / ly
    return (dy <= 0.5) & (dx <= 0.5)


def square_c(x, y, l):
    return rect_c(x, y, l, l)


def setup_grid(xsize, ysize, Dx, Dy):

    # Simulation size
    # xsize = 12.5
    # ysize = 12.5
    # Pixel size
    # Dx = 0.005
    # Dy = 0.005

    # Borders
    xl = xsize / 2.0
    yl = ysize / 2.0

    xs = np.arange(-xl, xl + Dx, Dx)
    ys = np.arange(-yl, yl + Dy, Dy)

    # print 'simulated images are', nxpoints, 'x', nypoints

    # grid
    xx, yy = np.meshgrid(xs, ys, sparse=False, indexing='xy')

    return xx, yy, xs, ys, xl, yl
