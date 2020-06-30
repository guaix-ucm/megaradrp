#
# Copyright 2019-2020 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import logging

import numpy

import megaradrp.instrument.focalplane as fp
from numina.frame.utils import copy_img


def subtract_sky(img, ignored_sky_bundles=None, logger=None):
    # Sky subtraction

    if logger is None:
        logger = logging.getLogger(__name__)

    logger.info('obtain fiber information')
    sky_img = copy_img(img)
    final_img = copy_img(img)
    fp_conf = fp.FocalPlaneConf.from_img(sky_img)
    # Sky fibers
    skyfibs = fp_conf.sky_fibers(valid_only=True,
                                   ignored_bundles=ignored_sky_bundles)
    logger.debug('sky fibers are: %s', skyfibs)
    # Create empty sky_data
    target_data = img[0].data

    target_map = img['WLMAP'].data
    sky_data = numpy.zeros_like(img[0].data)
    sky_map = numpy.zeros_like(img['WLMAP'].data)
    sky_img[0].data = sky_data

    for fibid in skyfibs:
        rowid = fibid - 1
        sky_data[rowid] = target_data[rowid]
        sky_map[rowid] = target_map[rowid]
    # Sum
    coldata = sky_data.sum(axis=0)
    colsum = sky_map.sum(axis=0)

    # Divide only where map is > 0
    mask = colsum > 0
    avg_sky = numpy.zeros_like(coldata)
    avg_sky[mask] = coldata[mask] / colsum[mask]

    # This should be done only on valid fibers
    logger.info('ignoring invalid fibers: %s', fp_conf.invalid_fibers())
    for fibid in fp_conf.valid_fibers():
        rowid = fibid - 1
        final_img[0].data[rowid, mask] = img[0].data[rowid, mask] - avg_sky[mask]
    # Update headers
    #
    return final_img, img, sky_img


def subtract_sky_rss(img, sky_img, ignored_sky_bundles=None, logger=None):
    """Subtract a sky image from an image"""
    # Sky subtraction

    if logger is None:
        logger = logging.getLogger(__name__)

    #logger.info('obtain fiber information')
    final_img = copy_img(img)
    # fiberconf_sky = dm.get_fiberconf(sky_img)
    # fiberconf_target = dm.get_fiberconf(img)

    logger.debug('using WLMAP extension to compute valid regions')

    v_map = img['WLMAP'].data > 0
    sky_map = numpy.zeros_like(img['WLMAP'].data)
    sky_data = sky_img[0].data
    sky_map[:] = v_map[:]

    # This should be done only on valid fibers
    #logger.info('ignoring invalid fibers: %s', fiberconf_target.invalid_fibers())
    final_img[0].data[v_map] = img[0].data[v_map] - sky_data[v_map]
    final_img[0].data[~v_map] = 0.0
    # Update headers
    #
    return final_img, img, sky_img