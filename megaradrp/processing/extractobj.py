
# Copyright 2011-2018 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Extract objects from RSS images"""

import numpy
import astropy.wcs
from scipy.spatial import KDTree


def extract_star(rssimage, position, npoints, fiberconf, logger=None):
    """Extract a star given its center and the number of fibers to extract"""

    # FIXME: handle several positions

    logger.info('extracting star')

    # fiberconf = datamodel.get_fiberconf(rssimage)
    logger.debug("Configuration UUID is %s", fiberconf.conf_id)

    rssdata = rssimage[0].data
    pdata = rssimage['wlmap'].data

    points = [position]
    fibers = fiberconf.conected_fibers(valid_only=True)
    grid_coords = []
    for fiber in fibers:
        grid_coords.append((fiber.x, fiber.y))
    # setup kdtree for searching
    kdtree = KDTree(grid_coords)

    # Other posibility is
    # query using radius instead
    # radius = 1.2
    # kdtree.query_ball_point(points, k=7, r=radius)

    dis_p, idx_p = kdtree.query(points, k=npoints)

    logger.info('Using %d nearest fibers', npoints)
    totals = []
    for diss, idxs, point in zip(dis_p, idx_p, points):
        # For each point
        logger.info('For point %s', point)
        colids = []
        coords = []
        for dis, idx in zip(diss, idxs):
            fiber = fibers[idx]
            colids.append(fiber.fibid - 1)
            logger.debug('adding fibid %d', fiber.fibid)

            coords.append((fiber.x, fiber.y))

        colids.sort()
        flux_fiber = rssdata[colids]
        flux_total = rssdata[colids].sum(axis=0)
        coverage_total = pdata[colids].sum(axis=0)

        max_cover = coverage_total.max()
        some_value_region = coverage_total > 0
        max_value_region = coverage_total == max_cover
        valid_region = max_value_region

        # Interval with maximum coverage
        nz_max, = numpy.nonzero(numpy.diff(max_value_region))
        # Interval with at least 1 fiber
        nz_some, = numpy.nonzero(numpy.diff(some_value_region))

        # Collapse the flux in the optimal region
        perf = flux_fiber[:, nz_max[0] + 1: nz_max[1] + 1].sum(axis=1)
        # Contribution of each fiber to the total flux, 1D
        perf_norm = perf / perf.sum()
        contributions = numpy.zeros(shape=(rssdata.shape[0],))
        contributions[colids] = perf_norm
        # Contribution of each fiber to the total flux, 2D
        flux_per_fiber = pdata * contributions[:, numpy.newaxis]
        flux_sum = flux_per_fiber.sum(axis=0)
        # In the region max_value_region, flux_sum == 1
        # In some_value_region 0 < flux_sum < 1
        # Outside is flux_sum == 0
        flux_correction = numpy.zeros_like(flux_sum)
        flux_correction[some_value_region] = 1.0 / flux_sum[some_value_region]

        # Limit to 10
        flux_correction = numpy.clip(flux_correction, 0, 10)

        flux_total_c = flux_total * flux_correction
        # plt.axvspan(nz_some[0], nz_some[1], alpha=0.2, facecolor='r')
        # plt.axvspan(nz_max[0], nz_max[1], alpha=0.2, facecolor='r')
        # plt.plot(flux_correction)
        # plt.show()
        #
        # plt.axvspan(nz_some[0], nz_some[1], alpha=0.2)
        # plt.axvspan(nz_max[0], nz_max[1], alpha=0.2)
        # plt.plot(flux_total)
        # plt.plot(flux_total * flux_correction, 'g')
        # ax2 = plt.gca().twinx()
        # ax2.plot(coverage_total)
        #
        # plt.show()
        pack = flux_total_c, colids, nz_max, nz_some
        # pack = flux_total_c
        totals.append(pack)

    return totals[0]


def compute_centroid(rssdata, fiberconf, c1, c2, point, logger=None):
    """Compute centroid near coordinates given by 'point'"""

    logger.debug("LCB configuration is %s", fiberconf.conf_id)

    fibers = fiberconf.conected_fibers(valid_only=True)
    grid_coords = []
    for fiber in fibers:
        grid_coords.append((fiber.x, fiber.y))
    # setup kdtree for searching
    kdtree = KDTree(grid_coords)

    # Other posibility is
    # query using radius instead
    # radius = 1.2
    # kdtree.query_ball_point(points, k=7, r=radius)

    npoints = 19
    # 1 + 6  for first ring
    # 1 + 6  + 12  for second ring
    # 1 + 6  + 12  + 18 for third ring
    points = [point]
    dis_p, idx_p = kdtree.query(points, k=npoints)

    logger.info('Using %d nearest fibers', npoints)
    for diss, idxs, point in zip(dis_p, idx_p, points):
        # For each point
        logger.info('For point %s', point)
        colids = []
        coords = []
        for dis, idx in zip(diss, idxs):
            fiber = fibers[idx]
            colids.append(fiber.fibid - 1)
            coords.append((fiber.x, fiber.y))

        coords = numpy.asarray(coords)
        flux_per_cell = rssdata[colids, c1:c2].mean(axis=1)
        flux_per_cell_total = flux_per_cell.sum()
        flux_per_cell_norm = flux_per_cell / flux_per_cell_total
        # centroid
        scf = coords.T * flux_per_cell_norm
        centroid = scf.sum(axis=1)
        logger.info('centroid: %s', centroid)
        # central coords
        c_coords = coords - centroid
        scf0 = scf - centroid[:, numpy.newaxis] * flux_per_cell_norm
        mc2 = numpy.dot(scf0, c_coords)
        logger.info('2nd order moments, x2=%f, y2=%f, xy=%f', mc2[0,0], mc2[1,1], mc2[0,1])
        return centroid


def compute_dar(img, datamodel, logger=None, debug_plot=False):
    """Compute Diferencial Atmospheric Refraction"""

    fiberconf = datamodel.get_fiberconf(img)
    wlcalib = astropy.wcs.WCS(img[0].header)

    rssdata = img[0].data
    cut1 = 500
    cut2 = 3500
    colids = []
    x = []
    y = []
    for fiber in fiberconf.fibers.values():
        colids.append(fiber.fibid - 1)
        x.append(fiber.x)
        y.append(fiber.y)

    cols = []
    xdar = []
    ydar = []
    delt = 50

    point = [2.0, 2.0]
    # Start in center of range
    ccenter = (cut2 + cut1) // 2
    # Left
    for c in range(ccenter, cut1, -delt):
        c1 = c - delt // 2
        c2 = c + delt // 2

        z = rssdata[colids, c1:c2].mean(axis=1)
        centroid = compute_centroid(rssdata, fiberconf, c1, c2, point, logger=logger)
        cols.append(c)
        xdar.append(centroid[0])
        ydar.append(centroid[1])
        point = centroid

    cols.reverse()
    xdar.reverse()
    ydar.reverse()

    point = [2.0, 2.0]
    # Star over
    # Right
    for c in range(ccenter, cut2, delt):
        c1 = c - delt // 2
        c2 = c + delt // 2
        z = rssdata[colids, c1:c2].mean(axis=1)
        centroid = compute_centroid(rssdata, fiberconf, c1, c2, point)
        cols.append(c)
        xdar.append(centroid[0])
        ydar.append(centroid[1])
        point = centroid

    rr = [[col, 0] for col in cols]
    world = wlcalib.wcs_pix2world(rr, 0)

    if debug_plot:
        import matplotlib.pyplot as plt
        import megaradrp.visualization as vis
        import megaradrp.simulation.refraction as r
        from astropy import units as u

        plt.subplots_adjust(hspace=0.5)
        plt.subplot(111)
        ax = plt.gca()
        plt.xlim([-8, 8])
        plt.ylim([-8, 8])
        col = vis.hexplot(ax, x, y, z, cmap=plt.cm.YlOrRd_r)
        plt.title("Fiber map, %s %s" % (c1, c2))
        cb = plt.colorbar(col)
        cb.set_label('counts')
        plt.show()

        zenith_distance = 60 * u.deg
        press = 79993.2 * u.Pa
        rel = 0.013333333
        temp = 11.5 * u.deg_C

        ll = r.differential_p(
                zenith_distance,
                wl=world[:,0] * u.AA,
                wl_reference=4025 * u.AA,
                temperature=temp,
                pressure=press,
                relative_humidity=rel,
        )
        plt.plot(world[:,0], xdar, '*-')
        plt.plot(world[:,0], ydar, '*-')
        plt.plot(world[:, 0], 2.0 * u.arcsec + ll.to(u.arcsec), '-')
        plt.show()

    return world[:, 0], xdar, ydar