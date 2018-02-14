#
# Copyright 2011-2017 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Base scientific recipe for MEGARA"""


import uuid
import math

from astropy.io import fits
import numpy as np
import numpy
import matplotlib.pyplot as plt

from numina.core import Parameter
from numina.core.requirements import ObservationResultRequirement
from numina.flow import SerialFlow
from numina.array import combine

from megaradrp.core.recipe import MegaraBaseRecipe
import megaradrp.requirements as reqs
from megaradrp.processing.combine import basic_processing_with_combination
from megaradrp.utils import copy_img
from megaradrp.processing.aperture import ApertureExtractor
from megaradrp.processing.wavecalibration import WavelengthCalibrator
from megaradrp.processing.fiberflat import Splitter, FlipLR, FiberFlatCorrector
from megaradrp.processing.twilight import TwilightCorrector
from megaradrp.processing.extractobj import compute_centroid, compute_dar


class ImageRecipe(MegaraBaseRecipe):
    """Base Image."""

    # Requirements  
    obresult = ObservationResultRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_bpm = reqs.MasterBPMRequirement()
    master_slitflat = reqs.MasterSlitFlatRequirement()
    master_wlcalib = reqs.WavelengthCalibrationRequirement()
    master_fiberflat = reqs.MasterFiberFlatRequirement()
    master_twilight = reqs.MasterTwilightRequirement()
    master_traces = reqs.MasterAperturesRequirement()
    extraction_offset = Parameter([0.0], 'Offset traces for extraction')
    ignored_sky_bundles = Parameter([], 'Ignore these sky bundles')
    master_sensitivity = reqs.SensitivityRequirement()
    reference_extinction = reqs.ReferenceExtinction()
    relative_threshold = Parameter(0.3, 'Threshold for peak detection')

    def base_run(self, rinput):

        # 2D reduction
        flow1 = self.init_filters(rinput, rinput.obresult.configuration)
        img = basic_processing_with_combination(rinput, flow1, method=combine.median)
        hdr = img[0].header
        self.set_base_headers(hdr)

        self.save_intermediate_img(img, 'reduced_image.fits')

        reduced2d = copy_img(img)

        # 1D, extraction, Wl calibration, Flat fielding
        reduced_rss = self.run_reduction_1d(img,
            rinput.master_traces, rinput.master_wlcalib,
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

        reduced_rss =  flow2(img)
        return reduced_rss

    def run_sky_subtraction(self, img, ignored_sky_bundles=None):
        # Sky subtraction
        self.logger.info('obtain fiber information')
        sky_img = copy_img(img)
        final_img = copy_img(img)
        fiberconf = self.datamodel.get_fiberconf(sky_img)
        # Sky fibers
        skyfibs = fiberconf.sky_fibers(valid_only=True,
                                       ignored_bundles=ignored_sky_bundles)
        self.logger.debug('sky fibers are: %s', skyfibs)
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
            if False:
                plt.plot(sky_data[rowid])
                plt.title("%d" % fibid)
                plt.show()
        # Sum
        coldata = sky_data.sum(axis=0)
        colsum = sky_map.sum(axis=0)

        # Divide only where map is > 0
        mask = colsum > 0
        avg_sky = numpy.zeros_like(coldata)
        avg_sky[mask] = coldata[mask] / colsum[mask]

        # This should be done only on valid fibers
        # The information of which fiber is valid
        # is in the tracemap, not in the header
        for fibid in fiberconf.valid_fibers():
            rowid = fibid - 1
            final_img[0].data[rowid, mask] = img[0].data[rowid, mask] - avg_sky[mask]

        return final_img, img, sky_img

    def read_wcs(self, hdr):
        crpix = hdr['CRPIX1']
        wlr0 = hdr['CRVAL1']
        delt = hdr['CDELT1']
        return crpix, wlr0, delt

    def compute_dar(self, img):
        import numpy.polynomial.polynomial as pol

        wl, xdar, ydar = compute_dar(img, self.datamodel, logger=self.logger)
        print('DAR, x:', pol.polyfit(wl, xdar, deg=3))
        print('DAR: y:', pol.polyfit(wl, ydar, deg=3))

    def centroid(self, rssdata, fiberconf, c1, c2, point):
        return compute_centroid(rssdata, fiberconf, c1, c2, point, logger=self.logger)

    def generate_sensitivity(self, final, spectrum, star_interp, extinc_interp, cover1, cover2, sigma=20.0):

        from astropy.wcs import WCS
        import matplotlib.pyplot as plt
        from scipy.ndimage.filters import uniform_filter, gaussian_filter

        crpix, wlr0, delt = self.read_wcs(final[0].header)

        wcsl = WCS(final[0].header)

        r1 = numpy.arange(final[0].shape[1])
        r2 = r1 * 0.0
        lm = numpy.array([r1, r2])
        wavelen_ = wcsl.all_pix2world(lm.T, 0.0)
        wavelen = wavelen_[:, 0]

        airmass = final[0].header['AIRMASS']
        exptime = final[0].header['EXPTIME']

        response_0 = spectrum / exptime
        valid = response_0 > 0
        # In magAB
        # f(Jy) = 3631 * 10^-0.4 mAB

        response_1 = 3631 * numpy.power(10.0, -0.4 * (star_interp(wavelen) + extinc_interp(wavelen) * airmass))
        r0max = response_0.max()
        r1max = response_1.max()
        r0 = response_0 / r0max
        r1 = response_1 / r1max

        pixm1, pixm2 = cover1
        pixr1, pixr2 = cover2

        max_valid = numpy.zeros_like(valid)
        max_valid[pixm1:pixm2 + 1] = True

        partial_valid = numpy.zeros_like(valid)
        partial_valid[pixr1:pixr2 + 1] = True

        valid = numpy.ones_like(response_0)
        valid[pixm2:] = 0
        valid[:pixm1+1] = 0

        sampling = int(2.0 / delt)
        # print('sampling is', sampling)

        # with numpy.errstate(invalid='ignore', divide='ignore'):
        #    ratio = r0 / r1

        pixf1, pixf2 = int(math.floor(pixm1 +  2* sigma)), int(math.ceil(pixm2 - 2 * sigma))
        # calc1 = [[pixm1, pixm2, pixr1, pixr2, pixf1, pixf2], [0, 0, 0, 0, 0, 0]]

        flux_valid = numpy.zeros_like(valid, dtype='bool')
        flux_valid[pixf1:pixf2 + 1] = True

        # lm2 = numpy.array(calc1)
        # wavelen_ = wcsl.all_pix2world(lm2.T, 0.0)

        r0_ens = gaussian_filter(r0, sigma=sigma)
        # mf = uniform_filter(r0, size=resolution)
        #
        # plt.plot(wavelen, r0, 'b*-')
        # plt.plot(wavelen, r0_ens, 'r*-')
        # plt.axvspan(lm3[0], lm3[1], facecolor='g', alpha=0.2)
        # plt.axvspan(lm3[2], lm3[3], facecolor='r', alpha=0.2)
        # plt.axvspan(lm3[4], lm3[5], facecolor='b', alpha=0.2)
        # #plt.plot(wavelen, r1)
        # plt.show()

        ratio2 = r0_ens / r1
        s_response = ratio2 * (r0max / r1max)

        # FIXME: include history
        sens = fits.PrimaryHDU(s_response, header=final[0].header)
        sens.header['uuid'] = str(uuid.uuid1())
        sens.header['tunit'] = ('Jy', "Final units")
        sens.header['PIXLIMF1'] = (pixf1 + 1, "Start of valid flux calibration")
        sens.header['PIXLIMF2'] = (pixf2 + 1, "End of valid flux calibration")
        sens.header['PIXLIMR1'] = (pixr1 + 1, 'Start of region with at least one fiber')
        sens.header['PIXLIMR2'] = (pixr2 + 1, 'End of region with at least one fiber')
        sens.header['PIXLIMM1'] = (pixm1 + 1, 'Start of region with all fibers')
        sens.header['PIXLIMM2'] = (pixm2 + 1, 'End of region with all fibers')
        return sens
