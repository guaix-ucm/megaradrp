#
# Copyright 2016-2021 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Focus Telescope Recipe for Megara"""

from __future__ import division, print_function

import numpy
from numina.array import combine
from numina.core.dataholders import Result, Parameter, Requirement
from numina.core.requirements import ObservationResultRequirement
from numina.exceptions import RecipeError

from megaradrp.instrument.focalplane import FocalPlaneConf
from megaradrp.recipes.scientific.base import ImageRecipe
import megaradrp.requirements as reqs
from megaradrp.processing.combine import basic_processing_with_combination_frames


class FocusTelescopeRecipe(ImageRecipe):
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
    master_apertures = reqs.MasterAperturesRequirement(alias='master_traces')
    extraction_offset = Parameter([0.0], 'Offset traces for extraction', accept_scalar=True)
    master_wlcalib = reqs.WavelengthCalibrationRequirement()
    position = Requirement(list, "Position of the reference object", default=(0, 0))
    # Products
    focus_table = Result(float)

    # @numina.core.validator.validate
    def run(self, rinput):
        # Basic processing
        self.logger.info('start focus telescope')

        obresult = rinput.obresult

        flow = self.init_filters(rinput, obresult.configuration)

        coors = rinput.position
        focus_t = 'M2UZ'

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
            raise RecipeError(f'We have only {len(image_groups)} different focus')

        all_images = {}
        for focus, frames in image_groups.items():
            self.logger.info('processing focus %s', focus)

            try:
                img = basic_processing_with_combination_frames(frames, flow, method=combine.median, errors=False)

                self.save_intermediate_img(img, f'focus2d-{focus}.fits')

                # 1D, extraction, Wl calibration, Flat fielding
                _, img1d = self.run_reduction_1d(
                    img,
                    rinput.master_apertures,
                    rinput.master_wlcalib,
                    rinput.master_fiberflat,
                    offset=rinput.extraction_offset
                )

                do_sky_subtraction = True
                if do_sky_subtraction:
                    self.logger.info('start sky subtraction')
                    final, origin, sky = self.run_sky_subtraction(
                        img1d, rinput.ignored_sky_bundles)
                    self.logger.info('end sky subtraction')
                else:
                    final = img1d
                    origin = final
                    sky = final

                self.save_intermediate_img(final, f'focus1d-{focus}.fits')

                self.logger.info('find lines and compute FWHM')
                star_rss_fwhm = self.run_on_image(final, coors)
                all_images[focus] = star_rss_fwhm

            except ValueError:
                self.logger.info('focus %s cannot be processed', focus)

        self.logger.info('fit FWHM of star')
        final = self.reorder_and_fit(all_images)
        self.logger.info('best focus is %s', final)

        self.logger.info('end focus telescope')
        return self.create_result(focus_table=final)

    def run_on_image(self, img, coors):
        """Extract spectra, find peaks and compute FWHM."""

        from scipy.spatial import KDTree

        fp_conf = FocalPlaneConf.from_img(img)
        self.logger.debug("LCB configuration is %s", fp_conf.conf_id)
        rssdata = img[0].data
        cut1 = 1000
        cut2 = 3000
        points = [(0, 0)] # Center of fiber 313
        fibers = fp_conf.connected_fibers(valid_only=True)
        grid_coords = []
        for fiber in fibers:
            grid_coords.append((fiber.x, fiber.y))
        # setup kdtree for searching
        kdtree = KDTree(grid_coords)

        # Other posibility is
        # query using radius instead
        # radius = 1.2
        # kdtree.query_ball_point(points, k=7, r=radius)

        npoints = 19 + 18
        # 1 + 6  for first ring
        # 1 + 6  + 12  for second ring
        # 1 + 6  + 12  + 18 for third ring
        dis_p, idx_p = kdtree.query(points, k=npoints)

        self.logger.info('Using %d nearest fibers', npoints)
        for diss, idxs, point in zip(dis_p, idx_p, points):
            # For each point
            self.logger.info('For point %s', point)
            colids = []
            coords = []
            for dis, idx in zip(diss, idxs):
                fiber = fibers[idx]
                colids.append(fiber.fibid - 1)
                coords.append((fiber.x, fiber.y))

            coords = numpy.asarray(coords)
            flux_per_cell = rssdata[colids, cut1:cut2].mean(axis=1)
            flux_per_cell_total = flux_per_cell.sum()
            flux_per_cell_norm = flux_per_cell / flux_per_cell_total
            # centroid
            scf = coords.T * flux_per_cell_norm
            centroid = scf.sum(axis=1)
            self.logger.info('centroid: %s', centroid)
            # central coords
            c_coords = coords - centroid
            scf0 = scf - centroid[:, numpy.newaxis] * flux_per_cell_norm
            mc2 = numpy.dot(scf0, c_coords)
            self.logger.info('2nd order moments, x2=%f, y2=%f, xy=%f', mc2[0,0], mc2[1,1], mc2[0,1])

        # FIXME: returning only 1 value for 1 star
        return mc2[0,0]

    def reorder_and_fit(self, all_images):
        """Fit all the values of FWHM to a 2nd degree polynomial and return minimum."""

        # We are assuming there is only 1 star
        focii = sorted(all_images.keys())

        ally = [all_images[focus] for focus in focii]

        try:
            res = numpy.polyfit(focii, ally, deg=2)
            self.logger.debug('fitting to deg 2 polynomial, done')
            self.logger.debug('parameters are %s', res)
            best = -res[1] / (2 * res[0])
        except ValueError as error:
            self.logger.warning("Error in fitting: %s", error)
            best = 0.0

        return best
