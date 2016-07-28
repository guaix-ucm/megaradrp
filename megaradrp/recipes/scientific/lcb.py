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
from base import ImageRecipe
import logging
import astropy.io.fits as fits
import numpy as np
from lcbmap import  Grid
from numina.core import Product
from numina.core.products import ArrayType
# from megaradrp.products import LCBCalibration
import re


_logger = logging.getLogger('numina.recipes.megara')


class LCBImageRecipe(ImageRecipe):
    """Process LCB images."""

    # master_lcb = Product(LCBCalibration)
    master_lcb = Product(ArrayType)

    def __init__(self):
        super(LCBImageRecipe, self).__init__()

    def run(self, rinput):

        _logger.info('starting LCB reduction')

        reduced, rssdata = super(LCBImageRecipe,self).run(rinput)

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
            _logger.info("vecinos: %s", vecinos)
            for elem in vecinos:
                neigh_info.append(grid.get_from_trace(elem))
                _logger.info('vecinos: %s', grid.get_from_trace(elem))

            centroid, sky, peak, sec_ord, cov = self.get_wcallib(100, 1000, vecinos,
                                                   rinput.wlcalib, rssdata,
                                                   np.array(neigh_info), grid)
            final_centroids.append(centroid)
            final_sky.append(sky)
            fiber.append(grid.get_fiber(centroid))
            peaks.append(peak)
            second_order.append(sec_ord)
            cova.append(cov)

        # master_lcb = self.generateJSON(points, final_centroids, final_sky, fiber, peaks, second_order, cova )

        tabla_final = self.generate_solution(points, final_centroids, final_sky, fiber, peaks, second_order, cova)

        return self.create_result(master_lcb=tabla_final)

    def get_wcallib(self, lambda1, lambda2, fibras, traces, rss, neigh_info, grid):

        # Take a look at == []
        indices = []
        wlcalib = []
        for elem in traces:
            if elem['aperture']['function']['coefficients']:
                wlcalib.append(elem['aperture']['function']['coefficients'])
                if len(indices)==0:
                    indices.append(0)
                else:
                    indices.append(indices[-1])
            else:
                indices.append(indices[-1]+1)

        wlcalib_aux = np.asarray(wlcalib)
        final, wcsdata = self.resample_rss_flux(rss, wlcalib_aux, indices)

        hdu_f = fits.PrimaryHDU(final)
        hdu_f.writeto('resample_rss.fits', clobber=True)

        fibras.sort()

        suma = np.sum(final[fibras,lambda1:lambda2],axis=1)
        sky = np.min(suma)
        sumaparcial = suma - sky

        neigh_info = neigh_info[np.argsort(neigh_info[:, 1])]

        centroid_x = np.multiply(sumaparcial,neigh_info[:,4])
        centroid_x = np.sum(centroid_x, axis=0)

        centroid_y = np.multiply(sumaparcial,neigh_info[:,5])
        centroid_y = np.sum(centroid_y, axis=0)

        sumatotal = np.sum(sumaparcial, axis=0)
        _logger.info( "sumatotal: %s", sumatotal)

        second_order = []
        aux = np.sum(np.multiply(suma,(neigh_info[:,4] - np.mean(neigh_info[:,4]))**2),axis=0)
        second_order.append(np.divide(aux ,np.sum(suma, axis=0)))
        _logger.info("Second order momentum X: %s", second_order[0])

        aux = np.sum(np.multiply(suma,(neigh_info[:,5] - np.mean(neigh_info[:,5]))**2),axis=0)
        second_order.append(np.divide(aux ,np.sum(suma, axis=0)))
        _logger.info("Second order momentum Y: %s", second_order[1])

        aux = np.multiply(neigh_info[:,5] - np.mean(neigh_info[:,5]),neigh_info[:,4] - np.mean(neigh_info[:,4]))
        aux = np.sum(np.multiply(aux,suma))
        cov = np.divide(aux ,np.sum(suma, axis=0))
        _logger.info("Cov X,Y: %s", cov)

        centroid_x = np.divide(centroid_x, sumatotal)
        _logger.info( "centroid_x: %s", centroid_x)

        centroid_y = np.divide(centroid_y, sumatotal)
        _logger.info("centroid_y: %s", centroid_y)

        centroid = [centroid_x, centroid_y]

        peak = np.sum(final[grid.get_fiber(centroid),lambda1:lambda2],axis=0)

        return centroid, sky, peak, second_order, cov

    def generateJSON(self, points, centroid, sky, fiber, peaks, second_order, cova):
        '''
        '''

        _logger.info('start JSON generation')

        result = []
        for cont, value in enumerate(points):
            obj = {
                'points': value,
                'centroid': centroid[cont],
                'sky':sky[cont],
                'fiber': fiber[cont],
                'peak': peaks[cont],
                'second_order': second_order[cont],
                'covariance': cova[cont]
            }
            result.append(obj)

        _logger.info('end JSON generation')

        return result

    def generate_solution(self, points, centroid, sky, fiber, peaks, second_order, cova):
        result = []
        for cont, value in enumerate(points):
            lista = (value[0], value[1], centroid[cont][0],centroid[cont][1], sky[cont], fiber[cont], peaks[cont], second_order[cont][0], second_order[cont][1], cova[cont])
            result.append(lista)
        return np.array(result, dtype=[('x_point','float'),('y_point','float'),('x_centroid','float'),('y_centroid','float'), ('sky','float'),('fiber','int'),('peak','float'),('x_second_order','float'), ('y_second_order','float'), ('covariance','float') ])