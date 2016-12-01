

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


if __name__ == '__main__':

    import matplotlib.pyplot as plt
    import megaradrp.processing.datamodel as dm
    import pkgutil
    from six import StringIO
    import astropy.io.fits as fits

    data = pkgutil.get_data('megaradrp', 'lcb_default_header.txt')

    default_hdr = StringIO(data)
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
