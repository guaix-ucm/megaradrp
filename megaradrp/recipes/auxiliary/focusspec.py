#
# Copyright 2016 Universidad Complutense de Madrid
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

"""Focus Spectrograph Recipe for Megara"""

from __future__ import division, print_function

import logging

import numpy
from scipy.spatial import cKDTree

from numina.core import Requirement, Parameter
from numina.core.dataholders import Product
from numina.core import DataFrameType
from numina.core.products import ArrayType
from numina.core.requirements import ObservationResultRequirement
from numina.core.products import LinesCatalog
import numina.array.fwhm as fmod
from numina.array.wavecal.statsummary import sigmaG
# FIXME: still using findpeaks1D functions
from numina.array.peaks.findpeaks1D import findPeaks_spectrum
from numina.array.peaks.findpeaks1D import refinePeaks_spectrum
from numina.flow.processing import BiasCorrector
from numina.flow import SerialFlow

from megaradrp.processing.trimover import OverscanCorrector, TrimImage
from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.products import TraceMap
from megaradrp.requirements import MasterBiasRequirement
from megaradrp.core.processing import apextract_tracemap
import astropy.io.fits as fits


_logger = logging.getLogger('numina.recipes.megara')


class FocusSpectrographRecipe(MegaraBaseRecipe):
    """Process Focus images and find best focus."""

    # Requirements
    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    tracemap = Requirement(TraceMap, 'Trace information of the Apertures')
    lines_catalog = Requirement(LinesCatalog, 'Catalog of lines')
    polynomial_degree = Parameter(2, 'Polynomial degree of arc calibration')
    # Products
    focus_table = Product(ArrayType)
    focus_image = Product(DataFrameType)

    def __init__(self):
        super(FocusSpectrographRecipe, self).__init__("0.1.0")

    def generate_image(self, final):
        from scipy.spatial import cKDTree

        voronoi_points = numpy.array(final[:, [0, 1]])
        x = numpy.arange(2048 * 2)
        y = numpy.arange(2056 * 2)

        test_points = numpy.transpose([numpy.tile(x, len(y)),numpy.repeat(y, len(x))])

        voronoi_kdtree = cKDTree(voronoi_points)

        test_point_dist, test_point_regions = voronoi_kdtree.query(test_points, k=1)
        final_image = test_point_regions.reshape((4112,4096)).astype('float64')
        final_image[:,:] = final[final_image[:,:].astype('int64'),2]
        return (final_image)

    def run(self, rinput):
        # Basic processing
        _logger.info('start focus spectrograph')
        flow = self.create_flow(rinput.master_bias)

        ever = []
        focci = []
        for frame in rinput.obresult.images:
            hdulist = frame.open()

            focus_val = hdulist[0].header['focus']
            focci.append(focus_val)

            _logger.info('processing frame %s', frame)
            hdulist = flow(hdulist)
            _logger.info('extract fibers')
            rssdata = apextract_tracemap(hdulist[0].data, rinput.tracemap)
            _logger.info('find lines and compute FWHM')
            lines_rss_fwhm = self.run_on_image(rssdata, rinput.tracemap)

            ever.append(lines_rss_fwhm)

        _logger.info('pair lines in images')
        line_fibers = self.filter_lines(ever)

        _logger.info('fit FWHM of lines')
        final = self.reorder_and_fit(line_fibers, focci)

        focus_median = numpy.median(final[:,2])
        _logger.info('median focus value is %5.2f', focus_median)

        _logger.info('generate focus image')
        image = self.generate_image(final)
        hdu = fits.PrimaryHDU(image)
        hdulist = fits.HDUList([hdu])

        _logger.info('end focus spectrograph')
        return self.create_result(focus_table=final, focus_image=hdulist)

    def get_focus(self, points, final):
        a = final[:, [0, 1]]

        b = points

        index = numpy.where(numpy.all(a==b,axis=1))
        return final[index[0],[2]]


    def create_flow(self, master_bias):

        biasmap = master_bias.open()[0].data

        flows = [OverscanCorrector(), TrimImage(), BiasCorrector(biasmap=biasmap)]
        basicflow = SerialFlow(flows)
        return basicflow

    def run_on_image(self, rssdata, tracemap):
        """Extract spectra, find peaks and compute FWHM."""

        # Extract the polynomials
        # FIXME: a little hackish
        pols = [numpy.poly1d(t['fitparms']) for t in tracemap]

        nwinwidth = 5
        times_sigma = 50.0
        lwidth = 20
        fpeaks = {}
        # FIXME: We are using here only one in 10 fibers
        for idx in range(0, len(rssdata), 10):
            # sampling every 10 fibers...
            row = rssdata[idx, :]
            # FIXME: this doesnt work if there are missing fibers
            # FIXME: cover is set, etc...
            the_pol = pols[idx]
            # find peaks
            threshold = numpy.median(row) + times_sigma * sigmaG(row)

            ipeaks_int = findPeaks_spectrum(row, nwinwidth, threshold)

            ipeaks_float = refinePeaks_spectrum(row, ipeaks_int, nwinwidth)
            fpeaks[idx] = []
            for peak, peak_f in zip(ipeaks_int, ipeaks_float):
                qslit = row[peak-lwidth:peak+lwidth]
                peak_val, fwhm = fmod.compute_fwhm_1d_simple(qslit, lwidth)
                peak_on_trace = the_pol(peak)
                fpeaks[idx].append((peak_f, peak_on_trace, fwhm))
            _logger.debug('found %d peaks in fiber %d', len(fpeaks[idx]), idx)
        return fpeaks

    def filter_lines(self, data):
        """Match lines between different images """

        ntotal = len(data)
        center = ntotal // 2
        base = data[center]
        _logger.debug('use image #%d as reference', center)
        line_fibers = {}
        maxdis = 2.0
        _logger.debug('matching lines up to %3.1f pixels', maxdis)
        for fiberid in base:

            ref = numpy.array(base[fiberid])
            outref = len(ref)
            savelines = {}
            line_fibers[fiberid] = savelines
            for i in range(outref):
                savelines[i] = {}
                savelines[i]['basedata'] = {}
                savelines[i]['centers'] = []
            kdtree = cKDTree(ref[:,:2])
            for i in range(ntotal):
                if i == center:
                    for j in range(outref):
                        savelines[j]['basedata']['coordinates'] = tuple(ref[j,:2])
                        savelines[j]['centers'].append(ref[j,2])
                    continue
                comp = numpy.array(data[i][fiberid])

                qdis, qidx = kdtree.query(comp[:,:2], distance_upper_bound=maxdis)
                for compidx, lidx in enumerate(qidx):
                    if lidx < outref:
                        savelines[lidx]['centers'].append(comp[compidx,2])

            remove_groups = []

            for ir in savelines:
                if len(savelines[ir]['centers']) != ntotal:
                    remove_groups.append(ir)

            for ir in remove_groups:
                _logger.debug('remove group of lines %d in fiber %d',ir, fiberid)
                del savelines[ir]

        return line_fibers

    def reorder_and_fit(self, line_fibers, focii):
        """Fit all the values of FWHM to a 2nd degree polynomial and return minimum."""

        l = 0
        for i in line_fibers:
            for _ in line_fibers[i]:
                l += 1

        _logger.debug('there are %d groups of lines to fit', l)
        ally = numpy.zeros((len(focii), l))
        final = numpy.zeros((l, 3))
        l = 0
        for i in line_fibers:
            for j in line_fibers[i]:
                ally[:,l] = line_fibers[i][j]['centers']
                final[l, :2] = line_fibers[i][j]['basedata']['coordinates']
                l += 1

        res = numpy.polyfit(focii, ally, deg=2)
        _logger.debug('fitting to deg 2 polynomial, done')
        best= -res[1] / (2*res[0])
        final[:,2] = best

        return final