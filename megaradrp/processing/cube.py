#
# Copyright 2017 Universidad Complutense de Madrid
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

import numpy as np
import math

from scipy import signal
from scipy.interpolate import RectBivariateSpline
from megaradrp.simulation.convolution import hex_c, square_c, setup_grid
from megaradrp.utils import copy_img
from megaradrp.datamodel import MegaraDataModel


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


def calc_matrix_from_fiberconf(fiberconf):
    # This should be in FIBERCONF...
    spos1_x = []
    spos1_y = []

    # FIBER in LOW LEFT corner is 614
    for fiber in fiberconf.conected_fibers():
        spos1_x.append(fiber.x)
        spos1_y.append(fiber.y)

    spos1_x = np.asarray(spos1_x)
    spos1_y = np.asarray(spos1_y)

    # FIBER in corner is 614
    # FIXME: workaround
    for fiber in fiberconf.conected_fibers():
        if fiber.fibid == 614:
            minx, miny = fiber.x, fiber.y
            if fiber.x < -6:
                # arcsec
                ascale = 0.536916 # 0.443 * 1.212
                print('scale is arcsec')
            else:
                # mm
                ascale = 0.443
                print('scale is mm')
            break
    else:
        raise ValueError('No fiber 614')

    refx, refy = minx / ascale, miny / ascale
    rpos1_x = (spos1_x - minx) / ascale
    rpos1_y = (spos1_y - miny) / ascale
    r0l_1 = np.array([rpos1_x, rpos1_y])
    return r0l_1, (refx, refy)


def calc_grid(scale=1.0):

    G_TYPE = 2
    ncol = 27
    nrow = 21
    r0l = calc_matrix(nrow, ncol, grid_type=G_TYPE)
    # r0l = R0 @ sl
    spos_x = scale * (r0l[0] - r0l[0].max() / 2)
    spos_y = scale * (r0l[1] - r0l[1].max() / 2)

    return spos_x, spos_y


def hexgrid_extremes(r0l, target_scale):
    # geometry
    h_hex = 0.5
    r_hex = 1 / M_SQRT3
    a_hex = 0.5 * h_hex * r_hex
    ha_hex = 6 * a_hex # detR0 == ha
    # compute extremes of hexgrid to rectangular grid
    # with pixel size 'scale'
    x0min, y0min = r0l.min(axis=1)
    x0max, y0max = r0l.max(axis=1)
    y1min = y0min - h_hex
    y1max = y0max + h_hex
    x1min = x0min - r_hex
    x1max = x0max + r_hex

    j1min = int(math.floor(x1min / target_scale + 0.5))
    i1min = int(math.floor(y1min / target_scale + 0.5))
    j1max = int(math.ceil(x1max / target_scale - 0.5))
    i1max = int(math.ceil(y1max / target_scale - 0.5))
    return (i1min, i1max), (j1min, j1max)


def create_cube(r0l, zval, target_scale=1.0):
    # geometry
    h_hex = 0.5
    r_hex = 1 / M_SQRT3
    a_hex = 0.5 * h_hex * r_hex
    ha_hex = 6 * a_hex # detR0 == ha

    R1 = target_scale * np.array([[1.0 ,0], [0,1]]) # Unit scale
    detR1 = np.linalg.det(R1)

    # compute extremes of hexgrid to rectangular grid
    # with pixel size 'scale'

    (i1min, i1max), (j1min, j1max) = hexgrid_extremes(r0l, target_scale)

    # Rectangular grid
    mk1 = np.arange(i1min, i1max+1)
    mk2 = np.arange(j1min, j1max+1)
    crow = len(mk1)
    ccol = len(mk2)
    # Result image
    # Third axis
    zval2 = my_atleast_2d(zval)
    print('zval shape is', zval2.shape)
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

    hex_kernel = hex_c(xx, yy, rad=r_hex, ang=0.0)
    square_kernel = square_c(xx, yy, target_scale)
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
        dk[s[1] - i1min, s[0] - j1min] = np.sum(we[:, np.newaxis] * zval2, axis=0)
    # scale with areas
    dk *= (target_scale**2 / ha_hex)
    return dk


def create_cube_from_rss(rss, target_scale_arcsec=1.0):

    HEX_SCALE = 0.536916 # Arcseconds

    target_scale = target_scale_arcsec / HEX_SCALE
    print('target scale is', target_scale)

    rss_data = rss[0].data
    # Operate on non-SKY fibers

    datamodel = MegaraDataModel()

    fiberconf = datamodel.get_fiberconf(rss)
    conected = fiberconf.conected_fibers()
    rows = [conf.fibid - 1 for conf in conected]
    #
    region = rss_data[rows, :]

    # FIXME: workaround
    # Get FUNIT keyword
    hdr = rss['FIBERS'].header
    funit = hdr.get('FUNIT', 'arcsec')
    pscale = hdr.get('PSCALE', 1.212)
    if funit == 'mm':
        print('Coordinates in mm')
        coord_scale = pscale
    else:
        print('Coordinates in arc sec')
        coord_scale = 1.0
    print('PSCALE', pscale)
    if False:
        import matplotlib.pyplot as plt
        import megaradrp.visualization as vi

        #
        spos_x = np.array([conf.x for conf in conected]) * coord_scale # Arcsec
        spos_y = np.array([conf.y for conf in conected]) * coord_scale # Arcsec
        zval = region[:, 2000]
        scale = HEX_SCALE # Arcsec
        plt.subplots_adjust(hspace=0.5)
        plt.subplot(111)
        ax = plt.gca()
        ll = 6.5
        plt.xlim([-ll, ll])
        plt.ylim([-ll, ll])

        plt.plot(spos_x, spos_y, 'r.')
        plt.show()

        plt.subplots_adjust(hspace=0.5)
        plt.subplot(111)
        ax = plt.gca()
        ll = 6.5
        plt.xlim([-ll, ll])
        plt.ylim([-ll, ll])
        col = vi.hexplot(ax, spos_x, spos_y, zval, scale=scale, cmap=plt.cm.YlOrRd_r)
        # plt.title("Fiber map")
        # cb = plt.colorbar(col)
        # cb.set_label('counts')
        plt.show()

    print('region.shape', region.shape)
    # plt.imshow(region, origin='lower', interpolation='bicubic')
    # plt.show()
    # The scale can be 'mm' or 'arcsec'  and it should come from the header
    r0l, (refx, refy) = calc_matrix_from_fiberconf(fiberconf)
    print('REF', refx, refy)

    (i1min, i1max), (j1min, j1max) = hexgrid_extremes(r0l, target_scale)
    cube_data = create_cube(r0l, region[:, :], target_scale)

    cube = copy_img(rss)
    # Move axis to put WL first
    # so that is last in FITS
    print('cube_data.shape', cube_data.shape)
    # plt.imshow(cube_data[:, :, 0], origin='lower', interpolation='bicubic')
    # plt.show()

    cube[0].data = np.moveaxis(cube_data, 2, 0)
    cube[0].data.astype('float32')
    # Merge headers
    merge_wcs(rss['FIBERS'].header, rss[0].header, out=cube[0].header)
    # Update values of WCS
    # CRPIX1, CRPIX2
    # CDELT1, CDELT2
    # minx, miny

    crpix_x = -refx / target_scale - j1min + 1 # FITS pixels are 1-based
    crpix_y = -refy / target_scale - i1min + 1
    print('crpix', crpix_y, crpix_x)
    # Map the center of original field
    #
    #
    cube[0].header['CRPIX1'] = crpix_x
    cube[0].header['CRPIX2'] = crpix_y
    cube[0].header['CDELT1'] = target_scale_arcsec / (3600.0)
    cube[0].header['CDELT2'] = target_scale_arcsec / (3600.0)
    # 2D from FIBERS
    # WL from PRIMARY
    # done
    return cube


def merge_wcs(hdr_sky, hdr_spec, out=None):

    if out is None:
        hdr = hdr_spec.copy()
    else:
        hdr = out

    # Extend header for third axis
    c_crpix = 'Pixel coordinate of reference point'
    c_cunit = 'Units of coordinate increment and value'
    hdr.set('CUNIT1', comment=c_cunit, after='CDELT1')
    hdr.set('CUNIT2', comment=c_cunit, after='CUNIT1')
    hdr.set('CUNIT3', comment=c_cunit, after='CUNIT2')
    hdr.set('CRPIX2', value=1, comment=c_crpix, after='CRPIX1')
    hdr.set('CRPIX3', value=1, comment=c_crpix, after='CRPIX2')
    hdr.set('CDELT3', after='CDELT2')
    hdr.set('CTYPE3', after='CTYPE2')
    hdr.set('CRVAL3', after='CRVAL2')
    c_pc = 'Coordinate transformation matrix element'
    hdr.set('PC1_1', value=1.0, comment=c_pc, after='CRVAL3')
    hdr.set('PC1_2', value=0.0, comment=c_pc, after='PC1_1')
    hdr.set('PC1_3', value=0.0, comment=c_pc, after='PC1_2')
    hdr.set('PC2_1', value=0.0, comment=c_pc, after='PC1_3')
    hdr.set('PC2_2', value=1.0, comment=c_pc, after='PC2_1')
    hdr.set('PC2_3', value=0.0, comment=c_pc, after='PC2_2')
    hdr.set('PC3_1', value=0.0, comment=c_pc, after='PC2_3')
    hdr.set('PC3_2', value=0.0, comment=c_pc, after='PC3_1')
    hdr.set('PC3_3', value=1.0, comment=c_pc, after='PC3_2')

    # Mapping, which keyword comes from each header
    mappings = [('CRPIX3', 'CRPIX1', 0),
                ('CDELT3', 'CDELT1', 0),
                ('CRVAL3', 'CRVAL1', 0),
                ('CTYPE3', 'CTYPE1', 0),
                ('CRPIX1', 'CRPIX1', 1),
                ('CDELT1', 'CDELT1', 1),
                ('CRVAL1', 'CRVAL1', 1),
                ('CTYPE1', 'CTYPE1', 1),
                ('CUNIT1', 'CUNIT1', 1),
                ('PC1_1', 'PC1_1', 1),
                ('PC1_2', 'PC1_2', 1),
                ('CRPIX2', 'CRPIX2', 1),
                ('CDELT2', 'CDELT2', 1),
                ('CRVAL2', 'CRVAL2', 1),
                ('CTYPE2', 'CTYPE2', 1),
                ('CUNIT2', 'CUNIT2', 1),
                ('PC2_1', 'PC2_1', 1),
                ('PC2_2', 'PC2_2', 1),
                ]

    idem_keys = [
        'LONPOLE',
    #    'LATPOLE',
        'RADESYS',
    #    'EQUINOX'
    ]
    for key in idem_keys:
        mp = (key, key, 1)
        mappings.append(mp)

    hdr_in = {}
    hdr_in[0] = hdr_spec
    hdr_in[1] = hdr_sky

    for dest, orig, idx in mappings:
        hdr_orig = hdr_in[idx]
        hdr[dest] = (hdr_orig[orig], hdr_orig.comments[orig])

    return hdr


def _simulate(seeing_fwhm=1.0, hex_scale=0.536916):
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


def _demo():
    import matplotlib.pyplot as plt

    seeing_fwhm = 1.1  # arcsec
    print('simulation')
    zval = _simulate(seeing_fwhm)
    print('done')
    HEX_SCALE = 0.536916

    _visualization(zval, scale=HEX_SCALE)

    print('zval shape is', zval.shape)
    G_TYPE = 2
    ncol = 27
    nrow = 21
    r0l = calc_matrix(nrow, ncol, grid_type=G_TYPE)
    # HEX_SCALE = 0.536916

    result = create_cube(r0l, zval, target_scale=0.5)
    print('result shape is', result.shape)
    plt.imshow(result[:,:,0], origin='lower', interpolation='bicubic')
    plt.show()


def _visualization(zval, scale=1.0):

    import matplotlib.pyplot as plt
    import megaradrp.visualization as vi

    spos_x, spos_y = calc_grid(scale=scale)
    plt.subplots_adjust(hspace=0.5)
    plt.subplot(111)
    ax = plt.gca()
    ll = 6.1
    plt.xlim([-ll, ll])
    plt.ylim([-ll, ll])
    col = vi.hexplot(ax, spos_x, spos_y, zval, scale=scale, cmap=plt.cm.YlOrRd_r)
    # plt.title("Fiber map")
    # cb = plt.colorbar(col)
    # cb.set_label('counts')
    plt.show()


def main(args=None):
    import argparse
    import astropy.io.fits as fits

    # parse command-line options
    parser = argparse.ArgumentParser(prog='convert_rss_cube')
    # positional parameters
    parser.add_argument("rss",
                        help="RSS file with fiber traces",
                        type=argparse.FileType('r'))
    parser.add_argument('-p', '--pixel-size', type=float, default=0.3,
                        metavar='PIXEL_SIZE',
                        help="Pixel size in arc seconds")
    parser.add_argument('-o', '--outfile', default='cube.fits',
                        help="Name of the output cube file")

    args = parser.parse_args(args=args)

    target_scale = args.pixel_size # Arcsec

    print('target scale is', target_scale, 'arcsec')

    with fits.open(args.rss) as rss:
        cube = create_cube_from_rss(rss, target_scale)

    cube.writeto(args.outfile, clobber=True)


if __name__ == '__main__':

    main()
