#
# Copyright 2017 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

from __future__ import print_function

import numpy as np
import math

from scipy import signal
from scipy.interpolate import RectBivariateSpline
from megaradrp.simulation.convolution import hex_c, square_c, setup_grid


M_SQRT3 = math.sqrt(3)


def my_atleast_2d(*arys):
    res = []
    for ary in arys:
        ary = np.asanyarray(ary)
        if len(ary.shape) == 0:
            result = ary.reshape(1, 1)
        elif len(ary.shape) == 1:
            result = ary[:, np.newaxis]
        else:
            result = ary
        res.append(result)
    if len(res) == 1:
        return res[0]
    else:
        return res


def calc_matrix(nrow, ncol, grid_type=2):

    R0 = np.array([[M_SQRT3 / 2,0], [-0.5,1]]) # Unit scale

    kcol = []
    krow = []
    for i in range(ncol):
        if grid_type == 1:
            s = (i + i % 2) // 2 # depending on type
        else:
            s = i // 2
        for j in range(nrow):
            kcol.append(i)
            krow.append(j+s)

    sl = np.array([kcol, krow]) # x y
    r0l = np.dot(R0, sl)
    # r0l = R0 @ sl

    return r0l


def calc_grid(scale=1.0):

    G_TYPE = 2
    ncol = 27
    nrow = 21
    r0l = calc_matrix(nrow, ncol, grid_type=G_TYPE)
    # r0l = R0 @ sl

    spos_x = scale * (r0l[0] - r0l[0].max() / 2)
    spos_y = scale * (r0l[1] - r0l[1].max() / 2)

    return spos_x, spos_y


def simulate(seeing_fwhm=1.0, hex_scale=0.536916):
    # simulation tools
    from megaradrp.simulation.atmosphere import generate_gaussian_profile
    from megaradrp.simulation.actions import simulate_point_like_profile

    r_hex = 1 / M_SQRT3
    FIBRAD_ANG = r_hex * hex_scale

    fibrad = FIBRAD_ANG  # arcsec
    seeing_profile = generate_gaussian_profile(seeing_fwhm)
    psf = None
    fraction_of_flux = simulate_point_like_profile(seeing_profile, psf, fibrad)

    spos_x, spos_y = calc_grid(scale=hex_scale)

    offpos0 = spos_x - 1.5
    offpos1 = spos_y + 1.6
    f_o_f = fraction_of_flux(offpos0, offpos1)
    b = 0.2
    zval = b + f_o_f
    zval = np.tile(zval[:, None], (1, 1000))
    zval = np.random.normal(zval, 0.01)
    return zval


def create_cube(zval, target_scale=1.0):
    # geometry
    h_hex = 0.5
    r_hex = 1 / M_SQRT3
    a_hex = 0.5 * h_hex * r_hex
    ha_hex = 6 * a_hex # detR0 == ha

    scale = target_scale
    R1 = scale * np.array([[1.0 ,0], [0,1]]) # Unit scale
    detR1 = np.linalg.det(R1)

    # In[7]:

    G_TYPE = 2
    ncol = 27
    nrow = 21
    r0l = calc_matrix(nrow, ncol, grid_type=G_TYPE)

    # compute extremes of hexgrid to rectangular grid
    # with pixel size 'scale'
    x0min, y0min = r0l.min(axis=1)
    x0max, y0max = r0l.max(axis=1)
    y1min = y0min - h_hex
    y1max = y0max + h_hex
    x1min = x0min - r_hex
    x1max = x0max + r_hex

    j1min = int(math.floor(x1min / scale + 0.5))
    i1min = int(math.floor(y1min / scale + 0.5))
    j1max = int(math.ceil(x1max / scale - 0.5))
    i1max = int(math.ceil(y1max / scale - 0.5))

    # Rectangular grid

    mk1 = np.arange(i1min, i1max+1)
    mk2 = np.arange(j1min, j1max+1)

    crow = len(mk1)
    ccol = len(mk2)
    # Result image
    # Third axis
    zval2 = my_atleast_2d(zval)
    # print('zval shape is', zval2.shape)
    # disp axis is last axis...
    dk = np.zeros((crow, ccol, zval2.shape[-1]))
    # print('result shape is ', dk.shape)
    # r1k = R1 @ sk
    sk = np.flipud(np.transpose([np.tile(mk1, len(mk2)), np.repeat(mk2, len(mk1))]).T)  # x y
    r1k = np.dot(R1, sk)

    # Compute convolution of hex and rect kernels
    Dx = 0.005
    Dy = 0.005
    xsize = ysize = 3.0
    xx, yy, xs, ys, xl, yl = setup_grid(xsize, ysize, Dx, Dy)

    hex_kernel = hex_c(xx, yy, rad=1.0 / M_SQRT3, ang=0.0)
    square_kernel = square_c(xx, yy, scale)
    convolved = signal.fftconvolve(hex_kernel, square_kernel, mode='same')
    kernel = convolved *(Dx *Dy)  / (detR1)
    rbs = RectBivariateSpline(xs, ys, kernel)
    # done

    # Loop to compute integrals...
    # This could be faster
    # zval could be 2D
    for s, r in zip(sk.T, r1k.T):
        allpos = -(r0l - r[:, np.newaxis])
        we = np.abs((rbs.ev(allpos[1], allpos[0])))
        we[we<0] = 0.0
        # k = np.sum(we * zval)
        # dk[s[1]-i1min, s[0]-j1min] = np.sum(we * zval)
        dk[s[1] - i1min, s[0] - j1min] = np.sum(we[:, np.newaxis] * zval2, axis=0)
    # scale with areas
    dk *= (scale**2 / ha_hex)
    return dk


if __name__ == '__main__':

    import matplotlib.pyplot as plt
    import megaradrp.visualization as vi
    import astropy.io.fits as fits

    seeing_fwhm = 1.1  # arcsec
    print('simulation')
    zval = simulate(seeing_fwhm)
    print('done')
    if False:
        HEX_SCALE = 0.536916
        spos_x, spos_y = calc_grid(scale=HEX_SCALE)
        plt.subplots_adjust(hspace=0.5)
        plt.subplot(111)
        ax = plt.gca()
        l = 6.1
        plt.xlim([-l, l])
        plt.ylim([-l, l])
        col = vi.hexplot(ax, spos_x, spos_y, zval, scale=HEX_SCALE, cmap=plt.cm.YlOrRd_r)
        #plt.title("Fiber map")
        cb = plt.colorbar(col)
        #cb.set_label('counts')
        plt.show()

    print('zval shape is', zval.shape)
    result = create_cube(zval, target_scale=0.5)
    print('result shape is', result.shape)
    plt.imshow(result[:,:,0], origin='lower', interpolation='bicubic')
    plt.show()
    fits.writeto('cube.fits', result.T, clobber=True)
