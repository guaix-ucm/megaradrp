#
# Copyright 2011-2016 Universidad Complutense de Madrid
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

"""Calibration Recipes for Megara"""


import math

from numina.core import Product
from numina.array.offrot import fit_offset_and_rotation
import numpy as np

from megaradrp.types import ProcessedRSS, ProcessedFrame
from megaradrp.recipes.scientific.base import ImageRecipe
from megaradrp.recipes.scientific.lcbmap import Grid
from megaradrp.types import LCBCalibration


class AcquireMOSRecipe(ImageRecipe):
    """Process MOS images."""

    #master_mos_json = Product(LCBCalibration)
    #master_mos = Product(ArrayType)

    reduced = Product(ProcessedFrame)
    #final = Product(ProcessedRSS)
    final = Product(ProcessedRSS)
    #sky = Product(ProcessedRSS)

    def run(self, rinput):

        self.logger.info('starting AC MOS reduction')

        reduced2d, reduced1d = super(AcquireMOSRecipe, self).base_run(rinput)
        # rssdata = rss_data[0].data

        self.save_intermediate_img(reduced2d, 'reduced2d.fits')
        self.save_intermediate_img(reduced1d, 'reduced1d.fits')

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

                colids = []
                coords = []
                for fiber in bundle.fibers.values():
                    colids.append(fiber.fibid - 1)
                    coords.append((fiber.x - central_fiber.x,
                                   fiber.y - central_fiber.y)
                                  )
                flux_per_cell = rssdata[colids, cut1:cut2].mean(axis=1)

                flux_per_cell_total = flux_per_cell.sum()
                #print(colids)
                #print(flux_per_cell / flux_per_cell_total)
                #print(coords)
                self.logger.debug('Compute centroid: %s', central_fiber.fibid)
                centroid = np.dot(np.asarray(coords).T, flux_per_cell) / flux_per_cell_total
                # Compute second order moments
                # I have this done somewhere
                print(key, central_coords, central_coords + centroid)
                p1.append(central_coords)
                q1.append(central_coords + centroid)

        offset, rot = fit_offset_and_rotation(np.array(p1), np.array(q1))
        angle = math.atan2(rot[1, 0], rot[0, 0])
        self.logger.info('offset is %s', offset)
        self.logger.info('rot matrix is %s', rot)
        self.logger.info('rot angle %s', angle)

        return self.create_result(reduced=reduced2d, final=final)
