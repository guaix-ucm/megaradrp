#
# Copyright 2016-2017 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# Megara DRP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Megara DRP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Megara DRP.  If not, see <http://www.gnu.org/licenses/>.
#

"""Focus Telescope Recipe for Megara"""

from __future__ import division, print_function

import numpy
from numina.array import combine
from numina.core.dataholders import Product
from numina.core.requirements import Requirement, ObservationResultRequirement
from numina.exceptions import RecipeError

from megaradrp.core.recipe import MegaraBaseRecipe
import megaradrp.requirements as reqs
from megaradrp.processing.combine import basic_processing_with_combination_frames
from megaradrp.processing.aperture import ApertureExtractor


class FocusTelescopeRecipe(MegaraBaseRecipe):
    """Process telescope focus images and find best focus.

    This recipe process a set of focus images obtained in
    **Focus Telescope** mode and returns an estimation
    of the telescope best focus.

    See Also
    --------
    megaradrp.recipes.auxiliary.focusspec.FocusSpectrographRecipe:
                recipe to measure the focus of the spectrograph

    Notes
    -----
    Images provided in `obresult` are grouped by the value of their
    FOCUST keyword. Groups of images are trimmed and corrected from overscan,
    bad pixel mask (if `master_bpm` is not None), bias and dark current
    (if `master_dark` is not None). Each group is then stacked using the median.

    The result of the combination is saved as an intermediate result, named
    'focus2d-#focus.fits', with #focus being the value of the focus
    of each group. Apertures are extracted in each combined image, and the
    resulting RSS file is saved as an intermediate result, named
    'focus1-#focus.fits'.

    For each image, the FWHM of the object at `position` is computed.

    Then, the FWHM  is fitted to a 2nd degree polynomial, and the focus
    corresponding to its minimum is obtained and returned in `focus_table`

    """

    # Requirements
    obresult = ObservationResultRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_bpm = reqs.MasterBPMRequirement()
    master_traces = reqs.MasterTraceMapRequirement()
    master_wlcalib = reqs.WavelengthCalibrationRequirement()
    position = Requirement(list, "Position of the reference object", default=(0, 0))
    # Products
    focus_table = Product(float)

    # @numina.core.validator.validate
    def run(self, rinput):
        # Basic processing
        self.logger.info('start focus telescope')

        obresult = rinput.obresult

        flow = self.init_filters(rinput, obresult.configuration)

        coors = rinput.position
        focus_t = 'FOCUST'

        image_groups = {}
        self.logger.info('group images by focus')

        for idx, frame in enumerate(obresult.frames):
            with frame.open() as img:
                focus_val = img[0].header[focus_t]
                if focus_val not in image_groups:
                    self.logger.debug('new focus %s', focus_val)
                    image_groups[focus_val] = []
                self.logger.debug('image %s in group %s', img, focus_val)
                image_groups[focus_val].append(frame)

        if len(image_groups) < 2:
            raise RecipeError('We have only {} different focus'.format(len(image_groups)))

        all_images = {}
        for focus, frames in image_groups.items():
            self.logger.info('processing focus %s', focus)

            try:
                img = basic_processing_with_combination_frames(frames, flow, method=combine.median, errors=False)
                calibrator_aper = ApertureExtractor(rinput.master_traces, self.datamodel)

                self.save_intermediate_img(img, 'focus2d-%s.fits' % (focus,))
                img1d = calibrator_aper(img)
                self.save_intermediate_img(img1d, 'focus1d-%s.fits' % (focus,))

                self.logger.info('find lines and compute FWHM')
                star_rss_fwhm = self.run_on_image(img1d, coors)
                all_images[focus] = star_rss_fwhm

            except ValueError:
                self.logger.info('focus %s cannot be processed', focus)

        self.logger.info('fit FWHM of star')
        final = self.reorder_and_fit(all_images)

        self.logger.info('end focus telescope')
        return self.create_result(focus_table=final)

    def run_on_image(self, img, coors):
        """Extract spectra, find peaks and compute FWHM."""

        # TODO
        return 0.0

    def reorder_and_fit(self, all_images):
        """Fit all the values of FWHM to a 2nd degree polynomial and return minimum."""

        # We are assuming there is only 1 star
        focii = sorted(all_images.keys())

        ally = [all_images[focus] for focus in focii]

        try:
            res = numpy.polyfit(focii, ally, deg=2)
            self.logger.debug('fitting to deg 2 polynomial, done')
            best = -res[1] / (2 * res[0])
        except ValueError as error:
            self.logger.warning("Error in fitting: %s", error)
            best = 0.0

        return best
