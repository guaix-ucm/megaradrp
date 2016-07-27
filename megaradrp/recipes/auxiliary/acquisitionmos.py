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
import logging

from numina.core import Product
from numina.core.products import ArrayType
import astropy.io.fits as fits
import numpy as np
import re

from megaradrp.recipes.scientific.base import ImageRecipe
from megaradrp.recipes.scientific.lcbmap import Grid
from megaradrp.products import LCBCalibration

_logger = logging.getLogger('numina.recipes.megara')


class AcquireMOSRecipe(ImageRecipe):
    """Process MOS images."""

    master_mos_json = Product(LCBCalibration)
    master_mos = Product(ArrayType)

    def __init__(self):
        super(AcquireMOSRecipe, self).__init__()

    def run(self, rinput):

        _logger.info('starting MOS reduction')

        reduced, rssdata = super(AcquireMOSRecipe, self).run(rinput)

        fits.writeto('rss.fits', rssdata, clobber=True)
        reduced.writeto('reduced.fits')

        bundles = []

        for key, value in reduced[1].header.items():
            bun = re.match(r"BUN(\d+)_T", key)
            if bun:
                if 'STAR' in value:
                    bundles.append(bun.group(1))

        _logger.info("bundles: %s", bundles)

        headers_list = np.array(reduced[1].header.items())
        bundles = [i.lstrip('0') for i in bundles]

        final_centroids = []
        final_sky = []

        fibers = []
        peaks = []
        second_order = []
        cova = []

        for elem in bundles:
            _logger.info('Bundle: %s', elem)
            spaxels_x = []
            spaxels_y = []
            spaxels_b = []
            fiber = []

            for aux in headers_list[np.where(headers_list[:, 1] == elem)]:
                fibb = re.match(r"FIB(\d+)_B", aux[0])
                if fibb:
                    _logger.info('FIB%s_X: %s', fibb.group(1), reduced[1].header[
                        'FIB' + fibb.group(1) + '_X'])
                    _logger.info('FIB%s_Y: %s', fibb.group(1), reduced[1].header[
                        'FIB' + fibb.group(1) + '_Y'])

                    spaxels_x.append(
                        reduced[1].header['FIB' + fibb.group(1) + '_X'])
                    spaxels_y.append(
                        reduced[1].header['FIB' + fibb.group(1) + '_Y'])
                    fiber.append(int(fibb.group(1)))
                    spaxels_b.append(elem)

            spaxels = np.array(zip(spaxels_x, spaxels_y, spaxels_b, fiber)).astype(float)
            grid = Grid(spaxels)
            spaxels_aux = np.array(zip(fiber, spaxels_b, spaxels_x, spaxels_y)).astype(float)
            vecinos = spaxels[:, 0]

            centroid, sky, peak, sec_ord, cov = self.get_wcallib(100, 1000,
                                                                 vecinos,
                                                                 rinput.wlcalib,
                                                                 rssdata,
                                                                 spaxels_aux, grid)

            final_centroids.append(centroid)
            final_sky.append(sky)
            fibers.append(spaxels[grid.get_fiber(centroid),:][-1])
            peaks.append(peak)
            second_order.append(sec_ord)
            cova.append(cov)

        #FIXME: points should be taken from the coordenates. Calculate which ones
        points = [[0,0], [1,1], [2,2]]
        tabla_final = self.generate_solution(points, final_centroids,
                                             final_sky, fibers, peaks,
                                             second_order, cova)

        master_mos_json = self.generateJSON(points, final_centroids, final_sky, fibers, peaks, second_order, cova )

        return self.create_result(master_mos=tabla_final, master_mos_json=master_mos_json)
