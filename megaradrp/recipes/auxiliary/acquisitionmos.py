#
# Copyright 2011-2019 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Acquisition with the Fiber MOS"""


import math

import numpy as np
from numina.array.offrot import fit_offset_and_rotation
from numina.core import Result, Parameter
from numina.core.qc import QC

from megaradrp.datamodel import TargetType
from megaradrp.recipes.scientific.base import ImageRecipe
from megaradrp.types import ProcessedRSS, ProcessedFrame
from megaradrp.utils import add_collapsed_mos_extension


class AcquireMOSRecipe(ImageRecipe):
    """Process Acquisition MOS images.

    This recipe processes a set of acquisition images
    obtained in **MOS Acquisition** mode and returns
    the offset and rotation required to center the
    fiducial objects in their reference positions.

    See Also
    --------
    megaradrp.recipes.auxiliary.acquisitionlcb.AcquireLCBRecipe
    numina.array.offrot: Kabsch algorithm for offset and rotation

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
    `master_traces` and resampled according to the wavelength calibration in
    `master_wlcalib`. Then is divided by the `master_fiberflat`.
    The resulting RSS is saved as an intermediate
    result named 'reduced_rss.fits'. This RSS is also returned in the field
    `reduced_rss` of the recipe result.

    The sky is subtracted by combining the the fibers marked as `SKY`
    in the fibers configuration. The RSS with sky subtracted is returned ini the
    field `final_rss` of the recipe result.

    Then, the centroid of each fiducial object, marked as `REFERENCE` in the fibers
    configuration, is computed. The offset and rotation needed to center
    each fiducial object in its bundle is computed and returned

    """

    # Requirements are defined in base class
    extraction_region = Parameter(
        [1000, 3000],
        description='Region used to compute a mean flux',
        nelem=2
    )

    reduced_image = Result(ProcessedFrame)
    reduced_rss = Result(ProcessedRSS)
    final_rss = Result(ProcessedRSS)
    offset = Result(list)
    rotang = Result(float)

    def run(self, rinput):

        self.logger.info('starting AC MOS reduction')

        reduced2d, reduced1d = super(AcquireMOSRecipe, self).base_run(rinput)

        do_sky_subtraction = True
        if do_sky_subtraction:
            self.logger.info('start sky subtraction')
            isb = rinput.ignored_sky_bundles
            if isb:
                self.logger.info('sky bundles ignored: %s', isb)
            final, origin, sky = self.run_sky_subtraction(reduced1d,
                                                          ignored_sky_bundles=isb)
            self.logger.info('end sky subtraction')
        else:
            final =  reduced1d
            origin = final
            sky = final

        fiberconf = self.datamodel.get_fiberconf(final)

        cut1, cut2 = rinput.extraction_region

        self.logger.debug("MOS configuration is %s", fiberconf.conf_id)
        rssdata = final[0].data
        scale, funit = self.datamodel.fiber_scale_unit(final, unit=True)
        self.logger.debug('unit is %s', funit)
        platescale = self.datamodel.PLATESCALE

        p1 = []
        q1 = []
        temp = []
        for key, bundle in fiberconf.bundles.items():
            if bundle.target_type == TargetType.REFERENCE:
                self.logger.debug("%s %s %s", key, bundle.target_name, bundle.target_type)
                sorted_fibers = [bundle.fibers[key] for key in sorted(bundle.fibers)]
                central_fiber = sorted_fibers[3] # Central fiber is number 4 in the list
                central_coords = [central_fiber.x * scale, central_fiber.y * scale]
                #central_fiber_pair_id
                # Central fiber is
                self.logger.debug('Center fiber is %d', central_fiber.fibid)
                self.logger.debug('Center fiber coordinates %f %f arcsec', 
                             central_fiber.x * scale, central_fiber.y * scale)

                colids = []
                coords = []
                for fiber in sorted_fibers:
                    colids.append(fiber.fibid - 1)
                    coords.append((fiber.x, fiber.y))
                self.logger.debug('nearest fibers')
                self.logger.debug('%s', [col + 1 for col in colids])
                coords = np.asarray(coords) * scale
                flux_per_cell = rssdata[colids, cut1:cut2].mean(axis=1)
                flux_per_cell_total = flux_per_cell.sum()
                flux_per_cell_norm = flux_per_cell / flux_per_cell_total
                # centroid
                scf = coords.T * flux_per_cell_norm
                centroid = scf.sum(axis=1)
                self.logger.info('centroid: %s arcsec', list(centroid))
                self.logger.info('centroid: %s mm', list(centroid / platescale))
                # central coords
                c_coords = coords - centroid
                scf0 = scf - centroid[:, np.newaxis] * flux_per_cell_norm
                mc2 = np.dot(scf0, c_coords)
                self.logger.info('2nd order moments, x2=%f, y2=%f, xy=%f arcsec^2', mc2[0, 0], mc2[1, 1], mc2[0, 1])

                p1.append(central_coords)
                q1.append(centroid)
                temp.append((bundle.id, central_fiber.fibid, central_fiber.x * scale, central_fiber.y * scale, centroid[0], centroid[1]))

        if self.intermediate_results:
            with open("centroids.txt", "w") as fd:
                for entry in temp:
                    fd.write("%s %s %6.3f %6.3f %6.3f %6.3f\n" % entry)

        self.logger.info('compute offset and rotation with %d points', len(p1))
        if len(p1) == 0:
            self.logger.warn('cant compute offset and rotation with 0 points')
            offset = [0.0, 0.0]
            angle = 0.0
            qc = QC.BAD
        else:
            offset, rot = fit_offset_and_rotation(np.array(p1), np.array(q1))
            angle = math.atan2(rot[1, 0], rot[0, 0])
            angle = angle / math.pi * 180.0
            qc = QC.GOOD
            self.logger.info('offset is %s', offset)
            self.logger.info('rot matrix is %s', rot)
            self.logger.info('rot angle %5.2f deg', angle)

        final = add_collapsed_mos_extension(final)
        origin = add_collapsed_mos_extension(origin)

        return self.create_result(
            reduced_image=reduced2d,
            reduced_rss=reduced1d,
            final_rss=final,
            offset=offset,
            rotang=angle,
            qc=qc
        )
