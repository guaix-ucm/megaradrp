#
# Copyright 2017-2021 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

from __future__ import print_function

import math
from packaging.version import parse, Version

import matplotlib
import matplotlib.cbook as cbook
import matplotlib.collections as mcoll
import matplotlib.colors as mcolors
import matplotlib.transforms as mtrans
import matplotlib.transforms as mtransforms
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch

import numpy as np

import megaradrp.instrument.constants as cons
import megaradrp.processing.fixrss as fixrss


M_SQRT3 = math.sqrt(3)
M_1_SQRT3 = 1 / M_SQRT3


def hexplot(axis, x, y, z, scale=1.0, extent=None,
            cmap=None, norm=None, vmin=None, vmax=None,
            alpha=None, linewidths=None, edgecolors='none',
            **kwargs):
    """
    Make a hexagonal grid plot.

    Returns
    =======
    object: matplotlib.collections.PolyCollection

    """

    # I have to add these due to changes in private and protected interfaces
    mpl_version = parse(matplotlib.__version__)
    if mpl_version >= Version("3.4"):
        axis._process_unit_info([("x", x), ("y", y)], kwargs, convert=False)
    else:
        axis._process_unit_info(xdata=x, ydata=y, kwargs=kwargs)

    x, y, z = cbook.delete_masked_points(x, y, z)

    x = np.array(x, float)
    y = np.array(y, float)

    sx = 2 * M_1_SQRT3 * scale * 0.99
    sy = scale * 0.99

    if extent is not None:
        xmin, xmax, ymin, ymax = extent
    else:
        xmin, xmax = (np.amin(x - sx), np.amax(x + sx)) if len(x) else (0, 1)
        ymin, ymax = (np.amin(y - sy), np.amax(y + sy)) if len(y) else (0, 1)

        # to avoid issues with singular data, expand the min/max pairs
        xmin, xmax = mtrans.nonsingular(xmin, xmax, expander=0.1)
        ymin, ymax = mtrans.nonsingular(ymin, ymax, expander=0.1)

    padding = 1.0e-9 * (xmax - xmin)
    xmin -= padding
    xmax += padding

    n = len(x)
    polygon = np.zeros((6, 2), float)

    mx = my = 0.99 * scale
    polygon[:, 0] = mx * np.array([-0.5 * M_1_SQRT3, 0.5 * M_1_SQRT3, 1.0 * M_1_SQRT3,
                                   0.5 * M_1_SQRT3, -0.5 * M_1_SQRT3, -1.0 * M_1_SQRT3])
    polygon[:, 1] = my * np.array([0.5, 0.5, 0.0, -0.5, -0.5, 0.0])

    offsets = np.zeros((n, 2), float)
    offsets[:, 0] = x
    offsets[:, 1] = y

    # I have to add these due to changes in private and protected interfaces
    if mpl_version >= Version("3.3"):
        collection = mcoll.PolyCollection(
            [polygon],
            edgecolors=edgecolors,
            linewidths=linewidths,
            offsets=offsets,
            transOffset=mtransforms.AffineDeltaTransform(axis.transData)
        )
    else:
        collection = mcoll.PolyCollection(
            [polygon],
            edgecolors=edgecolors,
            linewidths=linewidths,
            offsets=offsets,
            transOffset=mtransforms.IdentityTransform(),
            offset_position="data"
        )

    if isinstance(norm, mcolors.LogNorm):
        if (z == 0).any():
            # make sure we have not zeros
            z += 1

    if norm is not None:
        if norm.vmin is None and norm.vmax is None:
            norm.autoscale(z)

    if norm is not None and not isinstance(norm, mcolors.Normalize):
        msg = "'norm' must be an instance of 'mcolors.Normalize'"
        raise ValueError(msg)

    collection.set_array(z)
    collection.set_cmap(cmap)
    collection.set_norm(norm)
    collection.set_alpha(alpha)
    collection.update(kwargs)

    if vmin is not None or vmax is not None:
        collection.set_clim(vmin, vmax)
    else:
        collection.autoscale_None()

    corners = ((xmin, ymin), (xmax, ymax))
    axis.update_datalim(corners)
    axis.autoscale_view(tight=True)

    # add the collection last
    axis.add_collection(collection, autolim=False)
    return collection


def transform_to_pixel(wcs, coords, ref=0):
    # We need only the WL axis
    coords_m = [coords, np.zeros_like(coords)]
    r = wcs.wcs_world2pix(np.transpose(coords_m), ref)
    first_col = r[:, 0].astype('int')
    return first_col


def is_inside(shape, *coords):

    for coord in coords:
        for shape1d, crds in zip(shape, coord):
            if not is_inside_1d(shape1d, crds):
                return False
    return True


def is_inside_1d(shape1d, coords):
    for coord in coords:
        if coord < 0:
            return False
        if coord >= shape1d:
            return False
    return True


def extract_region(wcs_wl, region, rssdata):
    cut1, cut2 = transform_to_pixel(wcs_wl, region)

    if is_inside_1d(rssdata.shape[1], [cut1, cut2]):
        z = rssdata[:, cut1:cut2].mean(axis=1)
        return z
    else:
        raise ValueError('region to average is outside image')


def extract_column(wcs_wl, col, rssdata):
    cut, = transform_to_pixel(wcs_wl, [col])

    if is_inside_1d(rssdata.shape[1], [cut]):
        # cut1, cut2 = args.average_region
        z = rssdata[:, cut]
        return z
    else:
        raise ValueError('column is outside image')


def extract_zval(rssdata, wcs_wl, coordinate_type, column=None, average_region=None, continuum_region=None):

    if continuum_region is None and average_region is None:
        raise ValueError('continuum_region and average_region are None')

    if rssdata.ndim == 1:
        # data is 1d
        print('RSS is 1D')
        zval = rssdata
    elif rssdata.ndim == 2:
        print('world coordinates in', coordinate_type)
        if column is not None:
            print('use column', column)
            zval = extract_column(wcs_wl, column, rssdata)
        else:
            print('compute average in', average_region)
            zval = extract_region(wcs_wl, average_region, rssdata)
        # Subtract continuum
        if continuum_region:
            print('compute continuum in', continuum_region)
            z0 = extract_region(wcs_wl, continuum_region, rssdata)
            zval -= z0
    else:
        raise ValueError('image must be 1D or 2D')
    return zval


def main(argv=None):
    """Process command line"""

    import argparse
    import json

    import matplotlib.pyplot as plt
    import astropy.io.fits as fits
    from astropy.wcs import WCS
    from astropy.visualization import simple_norm
    import astropy.units as u

    import megaradrp.datamodel as dm
    from megaradrp.instrument.focalplane import FocalPlaneConf
    from megaradrp.processing.cube import create_cube_from_rss
    from megaradrp.processing.wcs import update_wcs_from_ipa, compute_pa_from_ipa

    # scale of the LCB grid in mm
    scale_lcb = cons.SPAXEL_SCALE.to(u.mm).value

    parser = argparse.ArgumentParser(description='Display MEGARA RSS images')
    parser.add_argument('--wcs-grid', '--display-wcs-grid', action='store_true',
                        help='Display WCS grid')
    parser.add_argument('--wcs-pa-from-header', action='store_true',
                        help="Use PA angle from PC keys", dest='pa_from_header')
    parser.add_argument('--fix-missing', action='store_true',
                        help="Interpolate missing fibers")
    parser.add_argument('--average-region', nargs=2, default=[1000, 3000],
                        type=int, help='Region of the RSS averaged on display')
    parser.add_argument('--extname', '-e', default='PRIMARY',
                        help='Name of the extension used')
    parser.add_argument('--column', '-c', type=int,
                        help='Column of the RSS on display')
    parser.add_argument('--continuum-region', nargs=2,
                        type=int,
                        help='Region of the RSS used for continuum subtraction')
    parser.add_argument('--coordinate-type', choices=['pixel', 'wcs'],
                        default='pixel',
                        help='Types of coordinates used')
    parser.add_argument('--colormap', type=plt.get_cmap,
                        help='Name of a valid matplotlib colormap')
    parser.add_argument('--plot-sky', '--display-sky', action='store_true',
                        help='Display SKY bundles')
    parser.add_argument('--display-fibid', '--plot-fibid', action='store_true',
                        help='Display fiber IDs of the spaxels')
    parser.add_argument('--plot-nominal-config', action='store_true',
                        help='Plot nominal configuration, do not use the header')
    parser.add_argument('--hide-values', action='store_true',
                        help='Do not show values out of range')
    parser.add_argument('--title', help='Title of the plot')
    parser.add_argument('--label', help='Legend of the colorbar')
    parser.add_argument('--hex-size', type=float,
                        help=f'Size of the hexagons (default is {scale_lcb})',
                        default=scale_lcb)
    parser.add_argument('--hex-rel-size', type=float,
                        help='Scale the size of hexagons by a factor',
                        default=1.0)
    parser.add_argument('--min-cut', type=float,
                        help='Inferior cut level')
    parser.add_argument('--max-cut', type=float,
                        help='Superior cut level')
    parser.add_argument('--percent', type=float,
                        help='Compute cuts using percentiles')
    parser.add_argument('--stretch',
                        choices=['linear', 'sqrt', 'power', 'log', 'asinh'],
                        default='linear',
                        help='Name of the strech method used for display'
                        )

    group_c = parser.add_argument_group('contouring')
    group_c.add_argument('--contour-pixel-size', type=float, default=0.4,
                         help="Pixel size in arc seconds for image reconstruction")
    group_c.add_argument('--contour-levels',
                         help="Contour levels")
    group_c.add_argument('--contour', action='store_true',
                         help="Draw contours")
    group_c.add_argument('--contour-image',
                         help="Image for computing contours")
    group_c.add_argument('--contour-image-column', type=int,
                         help='Column of image used for contouring')
    group_c.add_argument('--contour-image-save',
                         help='Save image used for contouring')
    group_c.add_argument('--contour-image-region', nargs=2, default=[1000, 3000],
                         type=int, help='Region of the image used for contouring')
    group_c.add_argument('--contour-is-density', action='store_true',
                         help='The data is a magnitude that does not require scaling')

    parser.add_argument('rss', metavar='RSS', nargs='+',
                        help='RSS images to process')

    args = parser.parse_args(args=argv)

    for fname in args.rss:
        with fits.open(fname) as img:
            if args.fix_missing:
                fibid = 623
                print(f'interpolate fiber {fibid}')
                img = fixrss.fix_missing_fiber(img, fibid)
            if args.plot_nominal_config:
                insmode = img['FIBERS'].header['INSMODE']
                fp_conf = dm.get_fiberconf_default(insmode)
            else:
                fp_conf = FocalPlaneConf.from_img(img)
            plot_mask = np.ones((fp_conf.nfibers,), dtype=bool)
            if not args.plot_sky:
                skyfibers = fp_conf.sky_fibers()
                skyfibers.sort()
                skyfibers_idx = [(fibid - 1) for fibid in skyfibers]
                plot_mask[skyfibers_idx] = False

            x0 = np.empty((fp_conf.nfibers,))
            y0 = np.empty((fp_conf.nfibers,))
            num = {}
            names = {}
            # Key is fibid
            for _, fiber in sorted(fp_conf.fibers.items()):
                idx = fiber.fibid - 1
                x0[idx] = fiber.x
                y0[idx] = fiber.y
                num[idx] = fiber.fibid
                names[idx] = fiber.name

            extname = args.extname
            if args.coordinate_type == 'wcs':
                # read WCS from extname
                wcs_wl = WCS(img[extname].header)
            else:
                wcs_wl = WCS(naxis=len(img[extname].shape))

            rssdata = np.squeeze(img[extname].data)

            zval0 = extract_zval(
                rssdata, wcs_wl, args.coordinate_type,
                args.column, args.average_region, args.continuum_region
            )
            fig = plt.figure()
            
            x = x0[plot_mask]
            y = y0[plot_mask]
            zval = zval0[plot_mask]
            projection = None

            if args.wcs_grid:
                print('display wcs grid')
                if not args.pa_from_header:
                    print('compute PA from IPA')
                    ipa = img[extname].header['IPA']
                    pa = compute_pa_from_ipa(ipa)
                    print('IPA is', ipa, ' PA is', pa)
                    hdr_fibers = update_wcs_from_ipa(img['FIBERS'].header, pa)
                else:
                    hdr_fibers = img['FIBERS'].header

                projection = WCS(hdr_fibers)

            # plt.subplots_adjust(hspace=0.5)
            ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=projection)

#            ax.set_xlim([-6.5, 6.5])
#            ax.set_ylim([-6, 6])

            if args.wcs_grid:
                ax.coords.grid(color='black', alpha=1.0, linestyle='solid')

            norm = simple_norm(zval, args.stretch,
                               min_cut=args.min_cut,
                               max_cut=args.max_cut,
                               percent=args.percent
                               )

            if args.hide_values:
                zdisp = zval.copy()
                if norm.vmin is not None:
                    zdisp[zval < norm.vmin] = np.nan
                if norm.vmax is not None:
                    zdisp[zval > norm.vmax] = np.nan
            else:
                zdisp = zval
            scale = args.hex_rel_size * args.hex_size
            col = hexplot(ax, x, y, zdisp, scale=scale, cmap=args.colormap, norm=norm)
            if args.display_fibid:
                # Plot spaxel labels
                for xx, yy, fid in zip(x, y, num):
                    tp1 = TextPath(
                        (xx-0.05, yy+0.16), f"{fid:03d}", size=0.045)
                    tp2 = TextPath(
                        (xx-0.05, yy-0.18), f"{names[fid]}", size=0.045)
                    ax.add_patch(PathPatch(tp1, zorder=2, lw=0, fc="black"))
                    ax.add_patch(PathPatch(tp2, zorder=2, lw=0, fc="black"))

            if args.title is not None:
                ax.set_title(args.title)

            cb = fig.colorbar(col)
            if args.label is not None:
                cb.set_label(args.label)

            if args.contour:
                target_scale_arcsec = args.contour_pixel_size

                if args.contour_image is not None:
                    with fits.open(args.contour_image) as hdu_con:
                        extname = args.extname

                        if args.coordinate_type == 'wcs':
                            # read WCS from extname
                            wcs_wl_c = WCS(hdu_con[extname].header)
                        else:
                            wcs_wl_c = WCS(naxis=len(hdu_con[extname].shape))

                        rssdata_con = np.squeeze(hdu_con[extname].data)

                        zval_c = extract_zval(rssdata_con,
                                              wcs_wl_c,
                                              args.coordinate_type,
                                              column=args.contour_image_column,
                                              average_region=args.contour_image_region,
                                              continuum_region=None)
                else:
                    zval_c = zval0

                # Build synthetic rss... for reconstruction
                primary = fits.PrimaryHDU(data=zval_c[:, np.newaxis], header=img[extname].header)
                synt = fits.HDUList([primary, img['FIBERS']])

                if args.contour_image_save is not None:
                    synt.writeto(args.contour_image_save)

                conserve_flux = not args.contour_is_density
                order = 1
                s_cube = create_cube_from_rss(synt, order, target_scale_arcsec, conserve_flux=conserve_flux)
                cube_wcs = WCS(s_cube[0].header).celestial
                px, py = cube_wcs.wcs.crpix
                interp = np.squeeze(s_cube[0].data)
                td = mtransforms.Affine2D().translate(-px, -py).scale(target_scale_arcsec, target_scale_arcsec)
                tt_d = td + ax.transData
                # im = ax.imshow(interp, alpha=0.9, cmap='jet', transform=tt_d)
                # im = ax.imshow(interp, alpha=0.9, cmap='jet', transform=ax.get_transform(cube_wcs))
                if args.contour_levels is not None:
                    levels = json.loads(args.contour_levels)
                    mm = ax.contour(interp, levels, transform=tt_d)
                else:
                    mm = ax.contour(interp, transform=tt_d)
                print('contour levels', mm.levels)

            plt.show()


if __name__ == '__main__':
    main()
