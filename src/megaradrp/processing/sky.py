#
# Copyright 2019-2026 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

import logging

import numpy

import megaradrp.instrument.focalplane as fp
from numina.frame.utils import copy_img


def _subtract_sky_data(img, sky_data, sky_valid, fibers):
    """Subtract the sky from the fibers of a RSS

    The sky is subtracted in the given fibers, where the fiber is
    covered (extension WLMAP) and the sky is valid; elsewhere in the
    coverage of the fibers the values are not changed. The pixels
    outside the coverage of the fibers are set to 0.

    Parameters
    ----------
    img : astropy.io.fits.HDUList
        RSS, with the extension WLMAP
    sky_data : numpy.ndarray
        Sky, a RSS array or a spectrum (subtracted from all the fibers)
    sky_valid : numpy.ndarray
        Pixels where the sky is valid, with the shape of `sky_data`
    fibers : list of int
        Fiber IDs of the fibers where the sky is subtracted

    Returns
    -------
    astropy.io.fits.HDUList
        A new RSS, with the sky subtracted
    """
    final_img = copy_img(img)
    covered = img["WLMAP"].data > 0
    rows = numpy.zeros(covered.shape[0], dtype=bool)
    rows[numpy.asarray(fibers, dtype=int) - 1] = True
    region = covered & sky_valid & rows[:, numpy.newaxis]
    sky_full = numpy.broadcast_to(sky_data, covered.shape)
    final_img[0].data[region] = img[0].data[region] - sky_full[region]
    final_img[0].data[~covered] = 0.0
    return final_img


def subtract_sky(img, ignored_sky_bundles=None, logger=None):
    """Subtract the average spectrum of the sky fibers of the RSS

    The sky is the average of the valid fibers of the SKY bundles, in
    each wavelength, of those covered in that wavelength. It is
    subtracted from the valid fibers (see `_subtract_sky_data`).

    Returns
    -------
    tuple
        The RSS with the sky subtracted, the input RSS and a RSS with
        the data of the sky fibers
    """

    if logger is None:
        logger = logging.getLogger(__name__)

    logger.info("obtain fiber information")
    sky_img = copy_img(img)
    fp_conf = fp.FocalPlaneConf.from_img(sky_img)
    # Sky fibers
    skyfibs = fp_conf.sky_fibers(valid_only=True, ignored_bundles=ignored_sky_bundles)
    logger.debug("sky fibers are: %s", skyfibs)
    # Create empty sky_data
    target_data = img[0].data

    target_map = img["WLMAP"].data
    sky_data = numpy.zeros_like(img[0].data)
    sky_map = numpy.zeros_like(img["WLMAP"].data)
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

    logger.info("ignoring invalid fibers: %s", fp_conf.invalid_fibers())
    final_img = _subtract_sky_data(img, avg_sky, mask, fp_conf.valid_fibers())
    return final_img, img, sky_img


def subtract_sky_rss(img, sky_img, ignored_sky_bundles=None, logger=None):
    """Subtract a sky RSS from a RSS, fiber by fiber

    The sky is valid where the sky RSS is covered (its extension WLMAP,
    if present). It is subtracted from the valid fibers (see
    `_subtract_sky_data`).

    Returns
    -------
    tuple
        The RSS with the sky subtracted, the input RSS and the sky RSS
    """

    if logger is None:
        logger = logging.getLogger(__name__)

    logger.debug("using WLMAP extension to compute valid regions")
    if "WLMAP" in sky_img:
        sky_valid = sky_img["WLMAP"].data > 0
    else:
        sky_valid = numpy.ones(sky_img[0].data.shape, dtype=bool)
    fp_conf = fp.FocalPlaneConf.from_img(img)
    logger.info("ignoring invalid fibers: %s", fp_conf.invalid_fibers())
    final_img = _subtract_sky_data(img, sky_img[0].data, sky_valid, fp_conf.valid_fibers())
    return final_img, img, sky_img
