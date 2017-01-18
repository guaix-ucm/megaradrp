

import numpy as np

import matplotlib.cbook as cbook

import matplotlib.collections as mcoll
import matplotlib.colors as mcolors
import matplotlib.transforms as mtransforms
import matplotlib.transforms as mtrans


def hexplot(axis, x, y, z, extent=None,
           cmap=None, norm=None, vmin=None, vmax=None,
           alpha=None, linewidths=None, edgecolors='none',
           **kwargs):

    if not axis._hold:
        axis.cla()

    axis._process_unit_info(xdata=x, ydata=y, kwargs=kwargs)

    x, y, z = cbook.delete_masked_points(x, y, z)

    x = np.array(x, float)
    y = np.array(y, float)

    # hardcoded
    sx = (2 * 0.465) * 0.99
    sy = (2 * 0.268) * 0.99

    if extent is not None:
        xmin, xmax, ymin, ymax = extent
    else:
        xmin, xmax = (np.amin(x-sx), np.amax(x+sx)) if len(x) else (0, 1)
        ymin, ymax = (np.amin(y-sy), np.amax(y+sy)) if len(y) else (0, 1)

        # to avoid issues with singular data, expand the min/max pairs
        xmin, xmax = mtrans.nonsingular(xmin, xmax, expander=0.1)
        ymin, ymax = mtrans.nonsingular(ymin, ymax, expander=0.1)

    padding = 1.e-9 * (xmax - xmin)
    xmin -= padding
    xmax += padding

    n = len(x)
    polygon = np.zeros((6, 2), float)
    polygon[:, 0] = sx * np.array([-0.5, 0.5, 1.0, 0.5, -0.5, -1.0]) / 3.0
    polygon[:, 1] = sy * np.array([0.5, 0.5, 0.0, -0.5, -0.5, 0.0])

    #S = math.sqrt(3) / 2
    #polygon[:, 0] = sx * np.array([-0.5, 0.5, 1.0, 0.5, -0.5, -1.0])
    #polygon[:, 1] = sy * np.array([S, S, 0.0, -S, -S, 0.0])

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


def _demo():
    import matplotlib.pyplot as plt
    import megaradrp.processing.datamodel as dm
    import pkgutil
    from six import StringIO
    import astropy.io.fits as fits

    data = pkgutil.get_data('megaradrp', 'lcb_default_header.txt')

    default_hdr = StringIO(data.decode(encoding='utf8'))
    hdr_fiber = fits.header.Header.fromfile(default_hdr)
    fiberconf = dm.read_fibers_extension(hdr_fiber)

    x = []
    y = []
    for fiber in fiberconf.fibers.values():
        x.append(-fiber.x)
        y.append(-fiber.y)

    z = np.random.normal(1000, 100, size=len(x))
    z[0:7] = 0
    plt.subplots_adjust(hspace=0.5)
    plt.subplot(111)
    ax = plt.gca()
    plt.xlim([-8, 8])
    plt.ylim([-8, 8])
    col = hexplot(ax, x, y, z, cmap=plt.cm.YlOrRd_r)
    plt.title("Fiber map")
    cb = plt.colorbar(col)
    cb.set_label('counts')
    plt.show()


if __name__ == '__main__':

    import sys
    import math

    import matplotlib.pyplot as plt
    import astropy.io.fits as fits
    from astropy.wcs import WCS

    import megaradrp.processing.datamodel as dm

    for fname in sys.argv[1:]:
        with fits.open(fname) as img:
            datamodel = dm.MegaraDataModel()
            fiberconf = datamodel.get_fiberconf(img)
            x = []
            y = []
            for fiber in fiberconf.fibers.values():
                x.append(fiber.x)
                y.append(fiber.y)
            cut1 = 1000
            cut2 = 3000
            z = img[0].data[:, cut1:cut2].mean(axis=1)

            wcs2 = WCS(img['FIBERS'].header)

            wcs3 = WCS(naxis=2)
            wcs3.wcs.crpix = [0, 0]
            wcs3.wcs.ctype = ['RA---TAN', 'DEC--TAN']
            wcs3.wcs.cdelt = [1.0 / 3600.0, 1.0/ 3600.0]
            wcs3.wcs.crval = [75.906999999999996,  74.848500000000001]
            # Rotation around (0,0)

            ang = 20.0 / 180 * math.pi
            cs = math.cos(ang)
            ss = math.sin(ang)
            wcs3.wcs.pc = [[cs, -ss], [ss, cs]]

            fig = plt.figure()

            #plt.subplots_adjust(hspace=0.5)
            ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcs3)
            ax.set_xlim([-6, 6])
            ax.set_ylim([-6, 6])

            ax.coords.grid(color='black', alpha=1.0, linestyle='solid')
            col = hexplot(ax, x, y, z, cmap=plt.cm.YlOrRd_r)

            plt.title("LCB")
            cb = plt.colorbar(col)
            cb.set_label('counts')
            plt.show()

