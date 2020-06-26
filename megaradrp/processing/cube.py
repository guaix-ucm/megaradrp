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
from scipy import signal
import astropy.units as u
import astropy.wcs
from numina.frame.utils import copy_img

from megaradrp.instrument.focalplane import FocalPlaneConf
# from megaradrp.datamodel import MegaraDataModel
from megaradrp.core.utils import atleast_2d_last
import megaradrp.processing.wcs as mwcs
import megaradrp.processing.hexspline as hspline
import megaradrp.instrument.constants as cons

# Normalized hexagon geometry
M_SQRT3 = math.sqrt(3)
H_HEX = 0.5
R_HEX = 1 / M_SQRT3
A_HEX = 0.5 * H_HEX * R_HEX
HA_HEX = 6 * A_HEX # detR0 == ha


# Size scale of the spaxel grid in arcseconds
HEX_SCALE = (cons.GTC_FC_A_PLATESCALE * cons.SPAXEL_SCALE).to(u.arcsec).value


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

    R0 = np.array([[M_SQRT3 / 2,0], [-0.5,1]]) # Unit scale

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

    sl = np.array([kcol, krow]) # x y
    r0l = np.dot(R0, sl)
    # r0l = R0 @ sl
    return r0l


def calc_matrix_from_fiberconf(fibersconf):
    """

    Parameters
    ----------
    fibersconf : megaradrp.instrument.focalplance.FocalPlabeConf

    Returns
    -------

    """

    # TODO: This should be in FIBERCONFS...
    spos1_x = []
    spos1_y = []
    for fiber in fibersconf.connected_fibers():
        spos1_x.append(fiber.x)
        spos1_y.append(fiber.y)
    spos1_x = np.asarray(spos1_x)
    spos1_y = np.asarray(spos1_y)

    # FIXME: workaround
    # FIBER in LOW LEFT corner is 614
    REFID = 614
    ref_fiber = fibersconf.fibers[REFID]
    minx, miny = ref_fiber.x, ref_fiber.y
    if ref_fiber.x < -6:
        # arcsec
        ascale = HEX_SCALE
        # print('fiber coordinates in arcsec')
    else:
        # mm
        ascale = cons.SPAXEL_SCALE.to(u.mm).value
        # print('fiber coordinates in mm')
    refx, refy = minx / ascale, miny / ascale
    rpos1_x = (spos1_x - minx) / ascale
    rpos1_y = (spos1_y - miny) / ascale
    r0l_1 = np.array([rpos1_x, rpos1_y])
    return r0l_1, (refx, refy)


def calc_grid(scale=1.0):
    """

    Parameters
    ----------
    scale : float

    Returns
    -------

    """

    G_TYPE = 2 # Values for MEGARA
    ncol = 27 #
    nrow = 21 #
    r0l = calc_matrix(nrow, ncol, grid_type=G_TYPE)
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


def create_cube(r0l, zval, p=1, target_scale=1.0):
    """

    Parameters
    ----------
    r0l
    zval
    p : {1, 2}
    target_scale : float, optional

    Returns
    -------

    Raises
    ------
    ValueError
        If `p` > 2

    """
    # geometry
    # Interpolation method. Allowed values are:
    # P = 1 NN
    # P = 2 Linear
    if p > 2:
        raise ValueError('p > 2 not implemented')

    R1 = target_scale * np.array([[1.0 ,0], [0,1]]) # Unit scale

    # compute extremes of hexgrid to rectangular grid
    # with pixel size 'scale'

    (i1min, i1max), (j1min, j1max) = hexgrid_extremes(r0l, target_scale)

    # Rectangular grid
    mk1 = np.arange(i1min, i1max + 1)
    mk2 = np.arange(j1min, j1max + 1)
    crow = len(mk1)
    ccol = len(mk2)
    # Result image
    # Add third last axis
    zval2 = atleast_2d_last(zval)
    # disp axis is last axis...
    dk = np.zeros((crow, ccol, zval2.shape[-1]))
    # print('result shape is ', dk.shape)
    # r1k = R1 @ sk
    sk = np.flipud(np.transpose([np.tile(mk1, len(mk2)), np.repeat(mk2, len(mk1))]).T)  # x y
    r1k = np.dot(R1, sk)

    # Prefiltering
    # For p = 1, prefilter coefficients with p = 1, coeff = 1
    # For p = 2, prefilter coefficients with p = 2, coeff = 1
    # No prefiltering in zval2 is required if p <= 2

    rbs = hspline.rescaling_kernel(p, scale=target_scale)

    # Loop to compute integrals...
    for s, r in zip(sk.T, r1k.T):
        allpos = -(r0l - r[:, np.newaxis])
        we = np.abs((rbs.ev(allpos[1], allpos[0])))
        dk[s[1] - i1min, s[0] - j1min] = np.sum(we[:, np.newaxis] * zval2, axis=0)

    # Postfiltering
    # For p = 1, final image in NN, postfilter coefficients with n = 1
    # For p = 2, final image is linear, postfilter coefficients with n = 3
    #
    if p == 1:
        # Coefficients post filtering to n = 2 * p - 1 == 1
        cpk = dk
        # Nearest-neighbor samples equal to coefficients
        img = cpk
    elif p == 2:
        # Coefficients post filtering to n = 2 * p - 1 == 3
        cpk = np.zeros_like(dk)
        # last axis
        for k in range(dk.shape[-1]):
            cpk[..., k] = signal.cspline2d(dk[..., k])
        # Linear samples equal to coefficients
        img = cpk
    else:
        raise ValueError('p > 2 not implemented')

    return img


def create_cube_from_array(rss_data, fiberconf, p=1, target_scale_arcsec=1.0, conserve_flux=True):
    """

    Parameters
    ----------
    rss_data
    fiberconf : megaradrp.instrument.focalplance.FocalPlaneConf
    p : {1, 2}
    target_scale_arcsec : float
    conserve_flux : bool

    Returns
    -------

    """

    target_scale = target_scale_arcsec / HEX_SCALE
    conected = fiberconf.connected_fibers()
    rows = [conf.fibid - 1 for conf in conected]

    rss_data = atleast_2d_last(rss_data)

    region = rss_data[rows, :]

    r0l, (refx, refy) = calc_matrix_from_fiberconf(fiberconf)
    cube_data = create_cube(r0l, region[:, :], p, target_scale)
    # scale with areas
    if conserve_flux:
        cube_data *= (target_scale ** 2 / HA_HEX)
    result = np.moveaxis(cube_data, 2, 0)
    result.astype('float32')
    return result


def create_cube_from_rss(rss, p=1, target_scale_arcsec=1.0, conserve_flux=True):
    """

    Parameters
    ----------
    rss
    p : {1, 2}
    target_scale_arcsec : float, optional
    conserve_flux : bool, optional

    Returns
    -------

    """

    target_scale = target_scale_arcsec / HEX_SCALE
    # print('target scale is', target_scale)

    rss_data = rss[0].data
    # Operate on non-SKY fibers

    fiberconf = FocalPlaneConf.from_img(rss)
    conected = fiberconf.connected_fibers()
    rows = [conf.fibid - 1 for conf in conected]
    #
    region = rss_data[rows, :]

    # FIXME: workaround
    # Get FUNIT keyword
    r0l, (refx, refy) = calc_matrix_from_fiberconf(fiberconf)

    (i1min, i1max), (j1min, j1max) = hexgrid_extremes(r0l, target_scale)
    cube_data = create_cube(r0l, region[:, :], p, target_scale)

    if conserve_flux:
        # scale with areas
        cube_data *= (target_scale ** 2 / HA_HEX)

    cube = copy_img(rss)
    # Move axis to put WL first
    # so that is last in FITS
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
    # After shifting the array
    # refpixel is -i1min, -j1min
    crpix_x = -refx / target_scale - j1min
    crpix_y = -refy / target_scale - i1min
    # Map the center of original field
    #
    #
    cube[0].header['CRPIX1'] = crpix_x
    cube[0].header['CRPIX2'] = crpix_y
    cube[0].header['CDELT1'] = -target_scale_arcsec / (3600.0)
    cube[0].header['CDELT2'] = target_scale_arcsec / (3600.0)
    # 2D from FIBERS
    # WL from PRIMARY
    # done
    return cube


def recompute_wcs(hdr):
    """Recompute the WCS rotations from IPA """
    ipa = hdr['IPA']
    pa = mwcs.compute_pa_from_ipa(ipa)
    print('IPA angle is:', ipa, 'PA angle is', math.fmod(pa, 360))
    x = hdr['PC1_1']
    y = hdr['PC1_2']
    print('PA from header is:', np.rad2deg(math.atan2(y, x)))
    return mwcs.update_wcs_from_ipa(hdr, pa)


def merge_wcs(hdr_sky, hdr_spec, out=None):
    """Merge sky WCS with spectral WCS"""
    if out is None:
        hdr = hdr_spec.copy()
    else:
        hdr = out

    allw = astropy.wcs.find_all_wcs(hdr_spec)
    for w in allw:
        ss = w.wcs.alt
        merge_wcs_alt(hdr_sky, hdr_spec, hdr, spec_suffix=ss)

    return hdr


def merge_wcs_alt(hdr_sky, hdr_spec, out, spec_suffix=''):
    """Merge sky WCS with spectral WCS"""

    hdr = out
    s = spec_suffix
    sf = s
    # Extend header for third axis
    c_crpix = 'Pixel coordinate of reference point'
    c_cunit = 'Units of coordinate increment and value'
    hdr.set('CUNIT1{}'.format(sf), comment=c_cunit, after='CDELT1{}'.format(sf))
    hdr.set('CUNIT2{}'.format(sf), comment=c_cunit, after='CUNIT1{}'.format(sf))
    hdr.set('CUNIT3{}'.format(sf), value='', comment=c_cunit, after='CUNIT2{}'.format(sf))
    hdr.set('CRPIX2{}'.format(sf), value=1, comment=c_crpix, after='CRPIX1{}'.format(sf))
    hdr.set('CRPIX3{}'.format(sf), value=1, comment=c_crpix, after='CRPIX2{}'.format(sf))
    hdr.set('CDELT3{}'.format(sf), after='CDELT2{}'.format(sf))
    hdr.set('CTYPE3{}'.format(sf), after='CTYPE2{}'.format(sf))
    hdr.set('CRVAL3{}'.format(sf), after='CRVAL2{}'.format(sf))
    c_pc = 'Coordinate transformation matrix element'
    hdr.set('PC1_1{}'.format(sf), value=1.0, comment=c_pc, after='CRVAL3{}'.format(sf))
    hdr.set('PC1_2{}'.format(sf), value=0.0, comment=c_pc, after='PC1_1{}'.format(sf))
    hdr.set('PC2_1{}'.format(sf), value=0.0, comment=c_pc, after='PC1_2{}'.format(sf))
    hdr.set('PC2_2{}'.format(sf), value=1.0, comment=c_pc, after='PC2_1{}'.format(sf))
    hdr.set('PC3_3{}'.format(sf), value=1.0, comment=c_pc, after='PC2_2{}'.format(sf))

    # Mapping, which keyword comes from each header
    mappings = [('CRPIX3', 'CRPIX1', s, 0),
                ('CDELT3', 'CDELT1', s, 0),
                ('CRVAL3', 'CRVAL1', s, 0),
                ('CTYPE3', 'CTYPE1', s, 0),
                ('CRPIX1', 'CRPIX1', '', 1),
                ('CDELT1', 'CDELT1', '', 1),
                ('CRVAL1', 'CRVAL1', '', 1),
                ('CTYPE1', 'CTYPE1', '', 1),
                ('CUNIT1', 'CUNIT1', '', 1),
                ('PC1_1', 'PC1_1', '', 1),
                ('PC1_2', 'PC1_2', '', 1),
                ('CRPIX2', 'CRPIX2', '', 1),
                ('CDELT2', 'CDELT2', '', 1),
                ('CRVAL2', 'CRVAL2', '', 1),
                ('CTYPE2', 'CTYPE2', '', 1),
                ('CUNIT2', 'CUNIT2', '', 1),
                ('PC2_1', 'PC2_1', '', 1),
                ('PC2_2', 'PC2_2', '', 1),
                ('LONPOLE', 'LONPOLE', '', 1),
                ('RADESYS', 'READESYS', '', 1),
                ('specsys', 'SPECSYS', s, 0),
                ('ssysobs', 'SSYSOBS', s, 0),
                ('velosys', 'VELOSYS', s, 0)
                ]

    hdr_in = {}
    hdr_in[0] = hdr_spec
    hdr_in[1] = hdr_sky

    for dest, orig, key, idx in mappings:
        hdr_orig = hdr_in[idx]
        korig = orig + key
        kdest = dest + sf
        try:
            hdr[kdest] = hdr_orig[korig], hdr_orig.comments[korig]
        except KeyError:
            # Ignoring errors. Copy only if keyword exists
            pass

    return hdr


def _simulate(seeing_fwhm=1.0, hex_scale=HEX_SCALE):
    # simulation tools
    from numina.instrument.simulation.atmosphere import generate_gaussian_profile
    from megaradrp.simulation.actions import simulate_point_like_profile

    FIBRAD_ANG = R_HEX * hex_scale

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

    _visualization(zval, scale=HEX_SCALE)

    print('zval shape is', zval.shape)
    G_TYPE = 2
    ncol = 27
    nrow = 21
    r0l = calc_matrix(nrow, ncol, grid_type=G_TYPE)

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

    methods = {'nn': 1, 'linear': 2}

    parser.add_argument("rss",
                        help="RSS file with fiber traces",
                        type=argparse.FileType('rb'))
    parser.add_argument('-p', '--pixel-size', type=float, default=0.3,
                        metavar='PIXEL_SIZE',
                        help="Pixel size in arc seconds")
    parser.add_argument('-o', '--outfile', default='cube.fits',
                        help="Name of the output cube file")
    parser.add_argument('-d', '--disable-scaling', action='store_true',
                        help="Disable flux conservation")
    parser.add_argument('-m', '--method', action='store', choices=['nn', 'linear'],
                        default='nn', help="Method of interpolation")
    parser.add_argument('--wcs-pa-from-header', action='store_true',
                        help="Use PA angle from header", dest='pa_from_header')

    args = parser.parse_args(args=args)

    target_scale = args.pixel_size # Arcsec
    p = methods[args.method]
    print('interpolation method is "{}"'.format(args.method))
    print('target scale is', target_scale, 'arcsec')
    conserve_flux = not args.disable_scaling

    with fits.open(args.rss) as rss:
        if not args.pa_from_header:
            # Doing it here so the change is propagated to
            # all alternative coordinates
            print('recompute WCS from IPA')
            rss['FIBERS'].header = recompute_wcs(rss['FIBERS'].header)
        cube = create_cube_from_rss(rss, p, target_scale, conserve_flux=conserve_flux)

    cube.writeto(args.outfile, overwrite=True)


if __name__ == '__main__':

    main()
