#
# Copyright 2017-2021 Universidad Complutense de Madrid
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
from scipy import signal
import astropy.units as u
import astropy.wcs
from numina.frame.utils import copy_img

from megaradrp.instrument.focalplane import FocalPlaneConf
# from megaradrp.datamodel import MegaraDataModel
from megaradrp.core.utils import atleast_2d_last
import megaradrp.processing.fixrss as fixrss
import megaradrp.processing.hexgrid as hg
import megaradrp.processing.hexspline as hspline
import megaradrp.instrument.constants as cons

# Helper function for equivalence conversion
GTC_PLATESCALE = u.plate_scale(cons.GTC_FC_A_PLATESCALE)
# Size scale of the spaxel grid in arcseconds
HEX_SCALE = cons.SPAXEL_SCALE.to(u.arcsec, GTC_PLATESCALE).value


def calc_matrix_from_fiberconf(fpconf, refid=614):
    """
    Compute hexagonal grid matrix from FocalPlaneConf

    Parameters
    ----------
    fpconf : megaradrp.instrument.focalplane.FocalPlaneConf
    refid : int
        fiber ID of reference fiber for grid coordinates
    Returns
    -------

    """

    # TODO: This should be in FIBERCONFS...
    spos1_x = []
    spos1_y = []
    for fiber in fpconf.connected_fibers():
        spos1_x.append(fiber.x)
        spos1_y.append(fiber.y)
    spos1_x = np.asarray(spos1_x)
    spos1_y = np.asarray(spos1_y)

    # FIBER in LOW LEFT corner is 614
    ref_fiber = fpconf.fibers[refid]
    minx, miny = ref_fiber.x, ref_fiber.y
    if fpconf.funit == 'arcsec':
        # arcsec
        ascale = HEX_SCALE
    else:
        # mm
        # fpconf.funit == 'mm'
        ascale = cons.SPAXEL_SCALE.to(u.mm).value
    ref = minx / ascale, miny / ascale
    rpos1_x = (spos1_x - minx) / ascale
    rpos1_y = (spos1_y - miny) / ascale
    r0l_1 = np.array([rpos1_x, rpos1_y])
    return r0l_1, ref


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

    rr1 = target_scale * np.array([[1.0, 0], [0, 1]])  # Unit scale

    # compute extremes of hexgrid to rectangular grid
    # with pixel size 'scale'

    (i1min, i1max), (j1min, j1max) = hg.hexgrid_extremes(r0l, target_scale)

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
    # r1k = rr1 @ sk
    sk = np.flipud(np.transpose([np.tile(mk1, len(mk2)), np.repeat(mk2, len(mk1))]).T)  # x y
    r1k = np.dot(rr1, sk)

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
    Create a cube array from a 2D or 1D array and focal plane configuration

    Parameters
    ----------
    rss_data
    fiberconf : megaradrp.instrument.focalplane.FocalPlaneConf
    p : {1, 2}
    target_scale_arcsec : float
    conserve_flux : bool

    Returns
    -------
    np.ndarray

    """

    target_scale = target_scale_arcsec / HEX_SCALE
    conected = fiberconf.connected_fibers()
    rows = [conf.fibid - 1 for conf in conected]

    rss_data = atleast_2d_last(rss_data)

    region = rss_data[rows, :]

    r0l, _ = calc_matrix_from_fiberconf(fiberconf)
    cube_data = create_cube(r0l, region[:, :], p, target_scale)

    if conserve_flux:
        # scale with areas
        cube_data *= (target_scale ** 2 / hg.HA_HEX)
    # Move axis to put WL first
    # so that is last in FITS
    result = np.moveaxis(cube_data, 2, 0)
    result.astype('float32')
    return result


def create_cube_from_rss(rss, p=1, target_scale_arcsec=1.0, conserve_flux=True):
    """
    Create a cube HDUlist from a RSS HDUList

    Parameters
    ----------
    rss : fits.HDUList
    p : {1, 2}
    target_scale_arcsec : float, optional
    conserve_flux : bool, optional

    Returns
    -------
    fits.HDUList
    """

    fiberconf = FocalPlaneConf.from_img(rss)
    result_arr = create_cube_from_array(
        rss[0].data, fiberconf, p=p,
        target_scale_arcsec=target_scale_arcsec,
        conserve_flux=conserve_flux
    )

    cube = copy_img(rss)
    cube[0].data = result_arr

    sky_header = rss['FIBERS'].header.copy()
    spec_header = rss[0].header
    # Update values of sky WCS
    # CRPIX1, CRPIX2
    # CDELT1, CDELT2
    # minx, miny
    # After shifting the array
    # refpixel is -i1min, -j1min
    target_scale = target_scale_arcsec / HEX_SCALE
    r0l, (refx, refy) = calc_matrix_from_fiberconf(fiberconf)
    (i1min, i1max), (j1min, j1max) = hg.hexgrid_extremes(r0l, target_scale)
    crpix_x = -refx / target_scale - j1min
    crpix_y = -refy / target_scale - i1min
    # Map the center of original field
    sky_header['CRPIX1'] = crpix_x
    sky_header['CRPIX2'] = crpix_y
    sky_header['CDELT1'] = -target_scale_arcsec / 3600.0
    sky_header['CDELT2'] = target_scale_arcsec / 3600.0

    # Merge headers
    # 2D from FIBERS
    # WL from PRIMARY
    merge_wcs(sky_header, spec_header, out=cube[0].header)

    # done
    return cube


def merge_wcs(hdr_sky, hdr_spec, out=None):
    """Merge sky WCS with spectral WCS

    Works only with main WCS and B WCS

    """
    if out is None:
        hdr = hdr_spec.copy()
    else:
        hdr = out

    allw = astropy.wcs.find_all_wcs(hdr_spec)
    wcsnames = [w.wcs.alt for w in allw]
    for ss in wcsnames:
        merge_wcs_alt(hdr_sky, hdr_spec, hdr, spec_suffix=ss)
    return hdr


def merge_wcs_alt(hdr_sky, hdr_spec, out, spec_suffix=' '):
    """Merge sky WCS with spectral WCS"""

    hdr = out
    if spec_suffix == ' ':
        sf = ''
    else:
        sf = spec_suffix
    # Extend header for third axis
    c_crpix = 'Pixel coordinate of reference point'
    c_cunit = 'Units of coordinate increment and value'
    wcsname_s = f'WCSNAME{sf}'
    if wcsname_s in hdr:
        prev = wcsname_s
    else:
        prev = f'CTYPE1{sf}'

    hdr.set(f'WCSAXES{sf}', value=3, before=prev)
    if sf != '':
        hdr.set(f'WCSNAME{sf}', value='', after='PC3_3')
        hdr.set(f'CTYPE1{sf}', value='', after=f'WCSNAME{sf}')
        hdr.set(f'CRPIX1{sf}', value=1.0, after=f'CTYPE1{sf}')
        hdr.set(f'CRVAL1{sf}', value=0.0, after=f'CRPIX1{sf}')
        hdr.set(f'CDELT1{sf}', value=1.0, after=f'CRVAL1{sf}')
    hdr.set(f'CUNIT1{sf}', value='deg', comment=c_cunit, after=f'CDELT1{sf}')
    hdr.set(f'CTYPE2{sf}', after=f'CUNIT1{sf}')
    if sf != '':
        hdr.set(f'CRPIX2{sf}', value=1.0, after=f'CTYPE2{sf}')
        hdr.set(f'CRVAL2{sf}', value=0.0, after=f'CRPIX2{sf}')
        hdr.set(f'CDELT2{sf}', value=1.0, after=f'CRVAL2{sf}')
    hdr.set(f'CUNIT2{sf}', value='deg', comment=c_cunit, after=f'CDELT2{sf}')
    hdr.set(f'CRPIX2{sf}', value=1, comment=c_crpix, after=f'CTYPE2{sf}')
    hdr.set(f'CTYPE3{sf}', after=f'CUNIT2{sf}')
    hdr.set(f'CRPIX3{sf}', value=1, comment=c_crpix, after=f'CTYPE3{sf}')
    hdr.set(f'CRVAL3{sf}', after=f'CRPIX3{sf}')
    hdr.set(f'CDELT3{sf}', after=f'CRVAL3{sf}')
    hdr.set(f'CUNIT3{sf}', comment=c_cunit, after=f'CDELT3{sf}')
    c_pc = 'Coordinate transformation matrix element'
    hdr.set(f'PC1_1{sf}', value=1.0, comment=c_pc, after=f'CUNIT3{sf}')
    hdr.set(f'PC1_2{sf}', value=0.0, comment=c_pc, after=f'PC1_1{sf}')
    hdr.set(f'PC2_1{sf}', value=0.0, comment=c_pc, after=f'PC1_2{sf}')
    hdr.set(f'PC2_2{sf}', value=1.0, comment=c_pc, after=f'PC2_1{sf}')
    hdr.set(f'PC3_3{sf}', value=1.0, comment=c_pc, after=f'PC2_2{sf}')

    # Mapping, which keyword comes from each header
    mappings = [('CRPIX3', 'CRPIX1', sf, 0),
                ('CDELT3', 'CDELT1', sf, 0),
                ('CRVAL3', 'CRVAL1', sf, 0),
                ('CTYPE3', 'CTYPE1', sf, 0),
                ('CRPIX1', 'CRPIX1', '', 1),
                ('CDELT1', 'CDELT1', '', 1),
                ('CRVAL1', 'CRVAL1', '', 1),
                ('CTYPE1', 'CTYPE1', '', 1),
                ('CUNIT3', 'CUNIT1', sf, 0),
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
                ('LATPOLE', 'LATPOLE', '', 1),
                ('RADESYS', 'RADESYS', '', 1),
                ('EQUINOX', 'EQUINOX', '', 1),
                ('WCSNAME', 'WCSNAME', sf, 0),
                ('specsys', 'SPECSYS', sf, 0),
                ('ssysobs', 'SSYSOBS', sf, 0),
                ('velosys', 'VELOSYS', sf, 0)
                ]

    hdr_in = dict()
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
    parser.add_argument('--fix-missing', action='store_true',
                        help="Interpolate missing fibers")

    args = parser.parse_args(args=args)

    target_scale = args.pixel_size  # Arcsec
    p = methods[args.method]
    print(f'interpolation method is "{args.method}"')
    print('target scale is', target_scale, 'arcsec')
    conserve_flux = not args.disable_scaling

    with fits.open(args.rss) as rss:
        if not args.pa_from_header:
            # Doing it here so the change is propagated to
            # all alternative coordinates
            print('recompute WCS from IPA')
            ipa = rss['PRIMARY'].header['IPA']
            rss['FIBERS'].header = fixrss.recompute_wcs(rss['FIBERS'].header, ipa=ipa)
        if args.fix_missing:
            fibid = 623
            print(f'interpolate fiber {fibid}')
            rss = fixrss.fix_missing_fiber(rss, fibid)

        cube = create_cube_from_rss(rss, p, target_scale, conserve_flux=conserve_flux)

    cube.writeto(args.outfile, overwrite=True)


if __name__ == '__main__':

    main()
