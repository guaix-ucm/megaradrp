#
# Copyright 2011-2017 Universidad Complutense de Madrid
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

"""Acquisition with the Fiber MOS"""


import math

from numina.core import Product
from numina.array.offrot import fit_offset_and_rotation
import numpy as np

from megaradrp.types import ProcessedRSS, ProcessedFrame
from megaradrp.recipes.scientific.base import ImageRecipe


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

    Then, the centroid of each fiducial object, marked as `STAR` in the fibers
    configuration, is computed. The offset and rotation needed to center
    each fiducial object in its bundle is computed and returned

    """

    # Requirements are defined in base class

    reduced_image = Product(ProcessedFrame)
    reduced_rss = Product(ProcessedRSS)
    final_rss = Product(ProcessedRSS)
    offset = Product(list)
    rotang = Product(float)

    def run(self, rinput):

        self.logger.info('starting AC MOS reduction')

        reduced2d, reduced1d = super(AcquireMOSRecipe, self).base_run(rinput)

        do_sky_subtraction = True
        if do_sky_subtraction:
            self.logger.info('start sky subtraction')
            final, origin, sky = self.run_sky_subtraction(reduced1d)
            self.logger.info('end sky subtraction')
        else:
            final =  reduced1d
            origin = final
            sky = final

        fiberconf = self.datamodel.get_fiberconf(final)

        cut1 = 1000
        cut2 = 3000
        self.logger.debug("MOS configuration is %s", fiberconf.conf_id)
        rssdata = final[0].data
        p1 = []
        q1 = []
        for key, bundle in fiberconf.bundles.items():
            if bundle.target_type == 'STAR':
                self.logger.debug("%s %s %s", key, bundle.target_name, bundle.target_type)
                mm = bundle.fibers.values()
                central_fiber = mm[3] # Central fiber is number 4 in the list
                central_coords = [central_fiber.x, central_fiber.y]
                #central_fiber_pair_id
                # Central fiber is
                self.logger.debug('Center fiber is %d', central_fiber.fibid)
                self.logger.debug('Center fiber coordinates %f %f', central_fiber.x, central_fiber.y)

                colids = []
                coords = []
                for fiber in bundle.fibers.values():
                    colids.append(fiber.fibid - 1)
                    coords.append((fiber.x, fiber.y))

                coords = np.asarray(coords)
                flux_per_cell = rssdata[colids, cut1:cut2].mean(axis=1)
                flux_per_cell_total = flux_per_cell.sum()
                flux_per_cell_norm = flux_per_cell / flux_per_cell_total
                # centroid
                scf = coords.T * flux_per_cell_norm
                centroid = scf.sum(axis=1)
                self.logger.info('centroid: %s', centroid)
                # central coords
                c_coords = coords - centroid
                scf0 = scf - centroid[:, np.newaxis] * flux_per_cell_norm
                mc2 = np.dot(scf0, c_coords)
                self.logger.info('2nd order moments, x2=%f, y2=%f, xy=%f', mc2[0, 0], mc2[1, 1], mc2[0, 1])

                p1.append(central_coords)
                q1.append(centroid)

        self.logger.info('compute offset and rotation with %d points', len(p1))
        offset, rot = fit_offset_and_rotation(np.array(p1), np.array(q1))
        angle = math.atan2(rot[1, 0], rot[0, 0])
        self.logger.info('offset is %s', offset)
        self.logger.info('rot matrix is %s', rot)
        self.logger.info('rot angle %s', angle)

        return self.create_result(
            reduced_image=reduced2d,
            reduced_rss=reduced1d,
            final_rss=final,
            offset=offset,
            rotang=angle,
        )
