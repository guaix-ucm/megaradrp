
import numpy as np
import math


FWHM_G = 2*math.sqrt(2*math.log(2))


def gauss_profile2(x, y, I, c, s):
    rr2 = (x - c[0])**2 + (y-c[1])**2
    intens = I*np.exp(-0.5 * rr2 / s**2)
    return intens


def gauss_profile2_base(x, y, s):
    rr2 = x**2 + y**2
    intens = np.exp(-0.5 * rr2 / s**2) / (2*math.pi*s**2)
    return intens


# Returns True in inside hexagon, rectangle and square
def hex_c(x, y, rad):
    d = rad * 2
    dx = np.abs(x) / d
    dy = np.abs(y) / d
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
    xx, yy = np.meshgrid(xs, ys, sparse=False)

    return xx, yy, xs, ys, xl, yl
