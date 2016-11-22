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

'''Calibration Recipes for Megara'''
from megaradrp.recipes.scientific.base import ImageRecipe
import logging
import astropy.io.fits as fits
import numpy as np
from megaradrp.recipes.scientific.lcbmap import Grid
from numina.core import Product
from numina.core.products import ArrayType
from megaradrp.types import LCBCalibration
import re


_logger = logging.getLogger('numina.recipes.megara')


class AcquireLCBRecipe(ImageRecipe):
    """Process LCB images."""

    master_lcb_json = Product(LCBCalibration)
    master_lcb = Product(ArrayType)

    def run(self, rinput):

        _logger.info('starting LCB reduction')

        reduced, rssdata = super(AcquireLCBRecipe,self).run(rinput)

        fits.writeto('rss.fits', rssdata, clobber=True)
        reduced.writeto('reduced.fits')

        spaxels_x = []
        spaxels_y = []
        spaxels_b = []
        fiber = []

        for key, value in reduced[1].header.items():
            if 'FIB' in key:
                if '_X' in key:
                    spaxels_x.append(value)
                elif '_Y' in key:
                    spaxels_y.append(value)
                    fiber.append(int(re.findall("[-+]?\d+[\.]?\d*", key)[0]))
                elif '_B' in key:
                    spaxels_b.append(int(value))

        spaxels = np.array(zip(spaxels_x,spaxels_y, spaxels_b, fiber))

        grid = Grid(spaxels)

        points = [[0,0]]
        final_centroids = []
        final_sky = []
        fiber = []
        peaks = []
        second_order = []
        cova = []
        for point in points:
            vecinos = grid.get_neighbours(point, radius=0.6)
            neigh_info = []
            _logger.info("neighbours: %s", vecinos)

            for elem in vecinos:
                neigh_info.append(grid.get_from_trace(elem))
                _logger.info('neighbours: %s', grid.get_from_trace(elem))

            centroid, sky, peak, sec_ord, cov = self.get_wcallib(100, 1000, vecinos,
                                                   rinput.wlcalib, rssdata,
                                                   np.array(neigh_info), grid)

            final_centroids.append(centroid)
            final_sky.append(sky)
            fiber.append(spaxels[grid.get_fiber(centroid),:][-1])
            peaks.append(peak)
            second_order.append(sec_ord)
            cova.append(cov)

        master_lcb_json = self.generateJSON(points, final_centroids, final_sky, fiber, peaks, second_order, cova )

        tabla_final = self.generate_solution(points, final_centroids, final_sky, fiber, peaks, second_order, cova)

        return self.create_result(master_lcb=tabla_final, master_lcb_json=master_lcb_json)