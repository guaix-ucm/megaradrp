#
# Copyright 2017-2018 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import math

import matplotlib.cbook as cbook
import matplotlib.collections as mcoll
import matplotlib.colors as mcolors
import matplotlib.transforms as mtrans
import matplotlib.transforms as mtransforms
import numpy as np

from megaradrp.processing.wcs import update_wcs_from_ipa, compute_pa_from_ipa


M_SQRT3 = math.sqrt(3)

PLATESCALE = 1.2120  # arcsec / mm
SCALE = 0.443  # mm from center to center, upwards


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
    if not axis._hold:
        axis.cla()

    axis._process_unit_info(xdata=x, ydata=y, kwargs=kwargs)

    x, y, z = cbook.delete_masked_points(x, y, z)

    x = np.array(x, float)
    y = np.array(y, float)

    sx = 2 / M_SQRT3 * scale * 0.99
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

    S = 1 / math.sqrt(3)
    mx = my = 0.99 * scale
    polygon[:, 0] = mx * np.array([-0.5 * S, 0.5 * S, 1.0 * S, 0.5 * S, -0.5 * S, -1.0 * S])
    polygon[:, 1] = my * np.array([0.5, 0.5, 0.0, -0.5, -0.5, 0.0])

    offsets = np.zeros((n, 2), float)
    offsets[:, 0] = x
    offsets[:, 1] = y

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


def transform_to_pixel(wcs, coords, ref=1):
    # We need only the WL axis
    coords_m = [coords, np.zeros_like(coords)]
    r = wcs.wcs_world2pix(np.transpose(coords_m), ref)
    first_col = r[:,0].astype('int')
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
    # print(args.average_region, cut1, cut2)

    if is_inside_1d(rssdata.shape[1], [cut1, cut2]):
        # cut1, cut2 = args.average_region
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


def main(argv=None):
    """Process command line"""

    import argparse

    import matplotlib.pyplot as plt
    import astropy.io.fits as fits
    from astropy.wcs import WCS
    from astropy.visualization import simple_norm

    import megaradrp.datamodel as dm

    parser = argparse.ArgumentParser(description='Display MEGARA RSS images')
    parser.add_argument('--wcs-grid', action='store_true',
                        help='Display WCS grid')
    parser.add_argument('--wcs-pa-from-header', action='store_true',
                        help="Use PA angle from PC keys", dest='pa_from_header')
    parser.add_argument('--average-region', nargs=2, default=[1000, 3000],
                        type=int)
    parser.add_argument('--extname', '-e', default='PRIMARY')
    parser.add_argument('--column', '-c', type=int)
    parser.add_argument('--continuum-region', nargs=2,
                        type=int)
    parser.add_argument('--coordinate-type', choices=['pixel', 'wcs'],
                        default='pixel')
    parser.add_argument('--colormap', type=plt.get_cmap)
    parser.add_argument('--title')
    parser.add_argument('--min-cut', type=float)
    parser.add_argument('--max-cut', type=float)
    parser.add_argument('--percent', type=float)
    parser.add_argument('--stretch',
                        choices=['linear', 'sqrt', 'power', 'log', 'asinh'],
                        default='linear'
                        )
    parser.add_argument('rss', metavar='RSS', nargs='+',
                        help='RSS images to process')

    args = parser.parse_args(args=argv)

    for fname in args.rss:
        with fits.open(fname) as img:
            extname = args.extname
            # image checks
            bunit = img[extname].header.get('BUNIT', 'counts')
            # insmode = img[extname].header.get('INSMODE', 'LCB')
            # Image shape

            if args.coordinate_type == 'wcs':
                # read WCS from extname
                wcs_wl = WCS(img[extname].header)
            else:
                wcs_wl = WCS(naxis=len(img[extname].shape))

            datamodel = dm.MegaraDataModel()
            fiberconf = datamodel.get_fiberconf(img)

            x = []
            y = []
            # Key is fibid
            for _, fiber in sorted(fiberconf.fibers.items()):
                x.append(fiber.x)
                y.append(fiber.y)

            rssdata = np.squeeze(img[extname].data)

            # Handle different dimensionality
            if rssdata.ndim == 1:
                # data is 1d
                print('RSS is 1D')
                zval = rssdata
            elif rssdata.ndim == 2:
                print('coordinates in', args.coordinate_type)
                if args.column is not None:
                    print('use column', args.column)
                    zval = extract_column(wcs_wl, args.column, rssdata)
                else:
                    print('compute average in', args.average_region)
                    zval = extract_region(wcs_wl, args.average_region,
                                          img[extname].data)
                # Subtract continuum
                if args.continuum_region:
                    print('compute continuum in', args.continuum_region)
                    z0 = extract_region(wcs_wl, args.continuum_region,
                                        img[extname].data)
                    zval -= z0
            else:
                raise ValueError('image must be 1D or 2D')

            fig = plt.figure()

            projection = None

            if args.wcs_grid:
                print('display wcs grid')
                if not args.pa_from_header:
                    print('compute PA from IPA')
                    ipa = img[extname].header['IPA']
                    pa = compute_pa_from_ipa(ipa)
                    print('IPA is', ipa, ' PA is', pa)
                    hdr = update_wcs_from_ipa(img['FIBERS'].header, pa)
                else:
                    hdr = img['FIBERS'].header

                projection = WCS(hdr)

            # plt.subplots_adjust(hspace=0.5)
            ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=projection)
            ax.set_xlim([-6.5, 6.5])
            ax.set_ylim([-6.5, 6.5])

            if args.wcs_grid:
                ax.coords.grid(color='black', alpha=1.0, linestyle='solid')

            norm = simple_norm(zval, args.stretch,
                               min_cut=args.min_cut,
                               max_cut=args.max_cut,
                               percent=args.percent
                               )

            col = hexplot(ax, x, y, zval, scale=SCALE, cmap=args.colormap, norm=norm)

            if args.title is not None:
                plt.title(args.title)

            cb = plt.colorbar(col)
            cb.set_label(bunit)

            plt.show()


if __name__ == '__main__':

    main()