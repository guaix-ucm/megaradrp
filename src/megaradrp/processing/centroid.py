#
# Copyright 2020 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Compute centroids in MEGARA RSS images """

import math
import logging

import numpy
from scipy.spatial import KDTree
from numina.constants import FWHM_G

from megaradrp.instrument.focalplane import FocalPlaneConf


_logger = logging.getLogger(__name__)


def find_brightest_spaxel(final, extraction_region):
    """Find the centroid around the brightest spaxel"""
    fp_conf = FocalPlaneConf.from_img(final)
    _logger.debug("LCB configuration is %s", fp_conf.conf_id)

    rssdata = final[0].data
    cut1, cut2 = extraction_region

    flux_per_cell_all = rssdata[:, cut1:cut2].mean(axis=1)

    max_cell = flux_per_cell_all.argmax() + 1
    max_fiber_ = fp_conf.fibers[max_cell]

    _logger.info("maximum flux in spaxel %d -- %s", max_cell, max_fiber_.name)
    return max_fiber_


def calc_centroid(final, extraction_region, point, nrings):
    """compute the centroid around point"""
    import megaradrp.instrument.constants as cons
    import megaradrp.datamodel as dm

    fp_conf = FocalPlaneConf.from_img(final)
    _logger.debug("LCB configuration is %s", fp_conf.conf_id)

    rssdata = final[0].data
    cut1, cut2 = extraction_region

    points = [point]

    fibers = fp_conf.connected_fibers(valid_only=True)

    grid_coords = []
    for fiber in fibers:
        grid_coords.append((fiber.x, fiber.y))
    # setup kdtree for searching
    kdtree = KDTree(grid_coords)

    # Other possibility is
    # query using radius instead
    # radius = 1.2
    # kdtree.query_ball_point(points, k=7, r=radius)
    _logger.debug('adding %d nrings', nrings)
    npoints = 1 + 3 * nrings * (nrings + 1)
    _logger.debug('adding %d fibers', npoints)

    dis_p, idx_p = kdtree.query(points, k=npoints)

    _logger.info('Using %d nearest fibers', npoints)

    scale, funit = dm.fiber_scale_unit(final, unit=True)
    _logger.debug('unit is %s', funit)

    platescale = cons.GTC_FC_A_PLATESCALE.value
    positions = []
    for diss, idxs, point in zip(dis_p, idx_p, points):
        # For each point
        value = [p * scale for p in point]
        value_mm = [(v / platescale) for v in value]
        _logger.info('For point %s arcsec', value)
        _logger.info('For point %s mm', value_mm)
        colids = []
        coords = []
        for dis, idx in zip(diss, idxs):
            fiber = fibers[idx]
            colids.append(fiber.fibid - 1)
            coords.append((fiber.x, fiber.y))
        _logger.debug('nearest fibers')
        _logger.debug('%s', [col + 1 for col in colids])
        coords = numpy.asarray(coords) * scale
        # flux_per_cell = flux_per_cell_all[colids]
        flux_per_cell = rssdata[colids, cut1:cut2].mean(axis=1)
        flux_per_cell_total = flux_per_cell.sum()
        flux_per_cell_norm = flux_per_cell / flux_per_cell_total
        # centroid
        scf = coords.T * flux_per_cell_norm
        centroid = scf.sum(axis=1)
        _logger.info('centroid: %s arcsec', list(centroid))
        _logger.info('centroid: %s mm', list(centroid / platescale))
        # central coords
        c_coords = coords - centroid
        scf0 = scf - centroid[:, numpy.newaxis] * flux_per_cell_norm
        mc2 = numpy.dot(scf0, c_coords)
        _logger.info('2nd order moments, x2=%f, y2=%f, xy=%f arcsec^2', mc2[0, 0], mc2[1, 1], mc2[0, 1])
        if (mc2[0, 0] > 0) and (mc2[1, 1] > 0):
            _logger.info('FWHM , x=%f, y=%f arcsec',
                             FWHM_G * math.sqrt(mc2[0, 0]),
                             FWHM_G * math.sqrt(mc2[1, 1])
                             )
        positions.append(centroid / platescale)
    return positions[0]


def calc_centroid_brightest(final, extraction_region, nrings):
    """Compute the centroid around the brightest spaxel"""
    max_fiber_ = find_brightest_spaxel(final, extraction_region)
    point = (max_fiber_.x, max_fiber_.y)
    # The brightest spaxel
    position = calc_centroid(final, extraction_region, point, nrings)
    return position