#
# Copyright 2011-2020 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Acquisition with LCB"""


import math

import numpy as np
from scipy.spatial import KDTree

from numina.core import Result, Parameter
from numina.core.validator import range_validator
from numina.constants import FWHM_G

from megaradrp.instrument.focalplane import FocalPlaneConf
from megaradrp.ntypes import ProcessedRSS, ProcessedImage
from megaradrp.processing.centroid import calc_centroid_brightest
from megaradrp.recipes.scientific.base import ImageRecipe


class AcquireLCBRecipe(ImageRecipe):
    """Process Acquisition LCB images.

    This recipe processes a set of acquisition images
    obtained in **LCB Acquisition** mode and returns
    the offset and rotation required to center the
    fiduciary object in its reference positions.

    See Also
    --------
    megaradrp.recipes.auxiliary.acquisitionmos.AcquireMOSRecipe

    Notes
    -----
    Images provided by `obresult` are trimmed and corrected
    from overscan, bad pixel mask (if `master_bpm` is not None),
    bias, dark current (if `master_dark` is not None) and
    slit-flat (if `master_slitflat` is not None).

    Images thus corrected are the stacked using the median.
    The result of the combination is saved as an intermediate result, named
    'reduced_image.fits'. This combined image is also returned in the field
    `reduced_image` of the recipe result.

    The apertures in the 2D image are extracted, using the information in
    `master_apertures` and resampled according to the wavelength calibration in
    `master_wlcalib`. Then is divided by the `master_fiberflat`.
    The resulting RSS is saved as an intermediate
    result named 'reduced_rss.fits'. This RSS is also returned in the field
    `reduced_rss` of the recipe result.

    The sky is subtracted by combining the the fibers marked as `SKY`
    in the fibers configuration. The RSS with sky subtracted is returned ini the
    field `final_rss` of the recipe result.

    Then, the centroid of the fiduciary object nearest to the center of the field
    is computed. The offset needed to center
    the fiduciary object in the center of the LCB is returned.

    """

    # Requirements are defined in base class
    points = Parameter([(0, 0)], "Coordinates")
    nrings = Parameter(3, "Number of rings to extract the star",
                       validator=range_validator(minval=1))
    extraction_region = Parameter(
        [1000, 3000],
        description='Region used to compute a mean flux',
        nelem=2
    )

    reduced_image = Result(ProcessedImage)
    reduced_rss = Result(ProcessedRSS)
    final_rss = Result(ProcessedRSS)
    offset = Result(list)
    rotang = Result(float)

    def run(self, rinput):

        self.logger.info('starting AC LCB reduction')

        reduced2d, reduced1d = super(AcquireLCBRecipe, self).base_run(rinput)
        # rssdata = rss_data[0].data

        do_sky_subtraction = True
        if do_sky_subtraction:
            self.logger.info('start sky subtraction')
            isb = rinput.ignored_sky_bundles
            if isb:
                self.logger.info('sky bundles ignored: %s', isb)
            final, origin, sky = self.run_sky_subtraction(
                reduced1d,
                sky_rss=rinput.sky_rss,
                ignored_sky_bundles=isb
            )
            self.logger.info('end sky subtraction')
        else:
            final =  reduced1d
            origin = final
            sky = final

        fp_conf = FocalPlaneConf.from_img(final)
        self.logger.debug("LCB configuration is %s", fp_conf.conf_id)

        centroid = calc_centroid_brightest(final, rinput.extraction_region, rinput.nrings)

        if False:
            self.compute_dar(final)

        return self.create_result(
            reduced_image=reduced2d,
            reduced_rss=reduced1d,
            final_rss=final,
            offset=-centroid
        )
