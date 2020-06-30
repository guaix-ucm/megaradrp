#
# Copyright 2011-2020 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Base scientific recipe for MEGARA"""


from numina.core import Parameter
from numina.core.requirements import ObservationResultRequirement
from numina.util.flow import SerialFlow
from numina.array import combine
from numina.frame.utils import copy_img

from megaradrp.core.recipe import MegaraBaseRecipe
import megaradrp.requirements as reqs
from megaradrp.processing.combine import basic_processing_with_combination

from megaradrp.processing.aperture import ApertureExtractor
from megaradrp.processing.wavecalibration import WavelengthCalibrator
from megaradrp.processing.fiberflat import FlipLR, FiberFlatCorrector
from megaradrp.processing.twilight import TwilightCorrector
from megaradrp.processing.extractobj import compute_centroid, compute_dar
from megaradrp.processing.sky import subtract_sky, subtract_sky_rss


class ImageRecipe(MegaraBaseRecipe):
    """Base Image."""

    # Requirements
    obresult = ObservationResultRequirement()
    method = Parameter(
        'median',
        description='Combination method',
        choices=['mean', 'median', 'sigmaclip']
    )
    method_kwargs = Parameter(
        dict(),
        description='Arguments for combination method',
        optional=True
    )
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_bpm = reqs.MasterBPMRequirement()
    master_slitflat = reqs.MasterSlitFlatRequirement()
    master_wlcalib = reqs.WavelengthCalibrationRequirement()
    master_fiberflat = reqs.MasterFiberFlatRequirement()
    master_twilight = reqs.MasterTwilightRequirement()
    master_apertures = reqs.MasterAperturesRequirement(alias='master_traces')
    sky_rss = reqs.SkyRSSRequirement(optional=True)
    extraction_offset = Parameter([0.0], 'Offset traces for extraction', accept_scalar=True)
    ignored_sky_bundles = Parameter([], 'Ignore these sky bundles')
    master_sensitivity = reqs.SensitivityRequirement()
    reference_extinction = reqs.ReferenceExtinction()
    relative_threshold = Parameter(0.3, 'Threshold for peak detection')
    diffuse_light_image = reqs.DiffuseLightRequirement()

    def base_run(self, rinput):

        # 2D reduction
        flow1 = self.init_filters(rinput, rinput.obresult.configuration)
        fmethod = getattr(combine, rinput.method)

        img = basic_processing_with_combination(
            rinput, flow1,
            method=fmethod,
            method_kwargs=rinput.method_kwargs
        )
        hdr = img[0].header
        self.set_base_headers(hdr)

        self.save_intermediate_img(img, 'reduced_image.fits')

        reduced2d = copy_img(img)

        # 1D, extraction, Wl calibration, Flat fielding
        reduced_rss = self.run_reduction_1d(img,
            rinput.master_apertures, rinput.master_wlcalib,
            rinput.master_fiberflat, rinput.master_twilight,
            offset=rinput.extraction_offset
        )
        self.save_intermediate_img(reduced_rss, 'reduced_rss.fits')

        return reduced2d, reduced_rss

    def run_reduction_1d(self, img, tracemap, wlcalib, fiberflat, twflat=None, offset=None):
        # 1D, extraction, Wl calibration, Flat fielding
        correctors = []
        correctors.append(ApertureExtractor(tracemap, self.datamodel, offset=offset))
        correctors.append(FlipLR())
        correctors.append(WavelengthCalibrator(wlcalib, self.datamodel))
        correctors.append(FiberFlatCorrector(fiberflat.open(), self.datamodel))

        if twflat:
            correctors.append(TwilightCorrector(twflat.open(), self.datamodel))

        flow2 = SerialFlow(correctors)

        reduced_rss = flow2(img)
        return reduced_rss

    def compute_dar(self, img):
        import numpy.polynomial.polynomial as pol

        wl, xdar, ydar = compute_dar(img, logger=self.logger)
        print('DAR, x:', pol.polyfit(wl, xdar, deg=3))
        print('DAR: y:', pol.polyfit(wl, ydar, deg=3))

    def centroid(self, rssdata, fiberconf, c1, c2, point):
        return compute_centroid(rssdata, fiberconf, c1, c2, point, logger=self.logger)

    def run_sky_subtraction(self, img, sky_rss=None, ignored_sky_bundles=None):

        if sky_rss is None:
            self.logger.info('compute sky from SKY bundles')
            if ignored_sky_bundles:
                self.logger.info('sky bundles ignored: %s', ignored_sky_bundles)
            return subtract_sky(img,
                                ignored_sky_bundles=ignored_sky_bundles,
                                logger=self.logger
                                )
        else:
            self.logger.info('use sky RSS image')
            return subtract_sky_rss(img, sky_img=sky_rss,
                                logger=self.logger
                                )