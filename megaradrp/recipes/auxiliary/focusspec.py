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
import numina.array.utils
import numina.array.fwhm as fmod
from numina.array.stats import robust_std as sigmaG
from numina.array.peaks.peakdet import find_peaks_indexes, refine_peaks
from numina.flow.processing import BiasCorrector
from numina.flow import SerialFlow
import astropy.io.fits as fits

from megaradrp.processing.trimover import OverscanCorrector, TrimImage
from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.products import WavelengthCalibration
from megaradrp.types import JSONstorage
import megaradrp.requirements as reqs
from megaradrp.core.processing import apextract_tracemap

_logger = logging.getLogger('numina.recipes.megara')

vph_thr = {'default':{'LR-I':{'min_distance':10,
                              'threshold':0.02},
                      'LR-R':{'min_distance':10,
                              'threshold':0.03},
                      'LR-V': {'min_distance':30,
                               'threshold':0.19},
                      'LR-Z': {'min_distance':60,
                               'threshold':0.02},
                      'LR-U':{'min_distance':10,
                              'threshold': 0.02,}
                      },
}


class FocusSpectrographRecipe(MegaraBaseRecipe):
    """Process Focus images and find best focus."""

    # Requirements
    obresult = ObservationResultRequirement()
    master_bias = reqs.MasterBiasRequirement()
    tracemap = reqs.MasterTraceMapRequirement()
    lines_catalog = Requirement(LinesCatalog, 'Catalog of lines')
    polynomial_degree = Parameter(2, 'Polynomial degree of arc calibration')
    wlcalib = Requirement(WavelengthCalibration,
                          'Wavelength calibration table')
    # Products
    focus_table = Product(ArrayType)
    focus_image = Product(DataFrameType)
    focus_wavelength = Product(JSONstorage)

    def __init__(self):
        super(FocusSpectrographRecipe, self).__init__("0.1.0")

    def generate_image(self, final):
        from scipy.spatial import cKDTree

        voronoi_points = numpy.array(final[:, [0, 1]])
        x = numpy.arange(2048 * 2)
        y = numpy.arange(2056 * 2)

        test_points = numpy.transpose(
            [numpy.tile(x, len(y)), numpy.repeat(y, len(x))])

        voronoi_kdtree = cKDTree(voronoi_points)

        test_point_dist, test_point_regions = voronoi_kdtree.query(test_points,
                                                                   k=1)
        final_image = test_point_regions.reshape((4112, 4096)).astype(
            'float64')
        final_image[:, :] = final[final_image[:, :].astype('int64'), 2]
        return (final_image)

    def run(self, rinput):
        # Basic processing
        _logger.info('start focus spectrograph')
        flow = self.create_flow(rinput.master_bias,
                                rinput.obresult.configuration.values)

        ever = []
        focci = []
        # FIXME: use tagger here
        vphs = []

        for frame in rinput.obresult.images:
            hdulist = frame.open()

            focus_val = hdulist[0].header['focus']
            vph_name = hdulist[0].header['VPH']
            focci.append(focus_val)

            _logger.info('processing frame %s', frame)
            hdulist = flow(hdulist)
            _logger.info('extract fibers')
            rssdata = apextract_tracemap(hdulist[0].data, rinput.tracemap)
            _logger.info('find lines and compute FWHM')
            lines_rss_fwhm = self.run_on_image(rssdata, rinput.tracemap, vph_name)

            ever.append(lines_rss_fwhm)

        _logger.info('pair lines in images')
        line_fibers = self.filter_lines(ever)

        focus_wavelength = self.generateJSON(ever,
                                             self.get_wlcalib(rinput.wlcalib),
                                             rinput.obresult.images)

        _logger.info('fit FWHM of lines')
        final = self.reorder_and_fit(line_fibers, focci)

        focus_median = numpy.median(final[:, 2])
        _logger.info('median focus value is %5.2f', focus_median)

        _logger.info('generate focus image')
        image = self.generate_image(final)
        hdu = fits.PrimaryHDU(image)
        hdulist = fits.HDUList([hdu])

        _logger.info('end focus spectrograph')
        return self.create_result(focus_table=final, focus_image=hdulist,
                                  focus_wavelength=focus_wavelength)

    def generateJSON(self, data, wlcalib, original_images):
        from numpy.polynomial.polynomial import polyval

        _logger.info('start JSON generation')

        result = {}
        counter = 0

        for image in data:
            name = original_images[counter].filename
            result[name] = {}
            for fiber, value in image.items():
                result[name][fiber] = []
                for arco in value:
                    try:
                        res = polyval(arco[0], wlcalib[fiber])
                        result[name][fiber].append(
                            [arco[0], arco[1], arco[2], res])
                    except:
                        _logger.error('Error in JSON generation. Check later...')
            counter += 1

        _logger.info('end JSON generation')

        return result

    def get_focus(self, points, final):
        a = final[:, [0, 1]]

        b = points

        index = numpy.where(numpy.all(a == b, axis=1))
        return final[index[0], [2]]

    def create_flow(self, master_bias, confFile):

        biasmap = master_bias.open()[0].data

        flows = [OverscanCorrector(confFile=confFile),
                 TrimImage(confFile=confFile),
                 BiasCorrector(biasmap=biasmap)]
        basicflow = SerialFlow(flows)
        return basicflow

    def run_on_image(self, rssdata, tracemap, current_vph):
        """Extract spectra, find peaks and compute FWHM."""

        # Extract the polynomials
        # FIXME: a little hackish
        valid_traces = get_valid_traces(tracemap)
        pols = [numpy.poly1d(t['fitparms']) for t in tracemap]

        vph_t = vph_thr['default']
        this_val = vph_t.get(current_vph)
        if this_val:
            flux_limit = this_val.get('flux_limit', 200000)
        else:
            flux_limit = 40000

        nwinwidth = 5
        times_sigma = 50.0
        lwidth = 20
        fpeaks = {}
        # FIXME: We are using here only one in 10 fibers
        for fibid in valid_traces[::10]:
            # sampling every 10 fibers...
            idx = fibid - 1
            row = rssdata[idx, :]

            the_pol = pols[idx]

            # FIXME: using here a different peak routine than in arc
            # find peaks
            threshold = numpy.median(row) + times_sigma * sigmaG(row)

            ipeaks_int1 = find_peaks_indexes(row, nwinwidth, threshold)
            # filter by flux
            _logger.info('Filtering peaks over %5.0f', flux_limit)
            ipeaks_vals = row[ipeaks_int1]
            mask = ipeaks_vals < flux_limit
            ipeaks_int = ipeaks_int1[mask]
            _logger.debug('LEN (ipeaks_int): %s', len(ipeaks_int))
            _logger.debug('ipeaks_int: %s', ipeaks_int)
            ipeaks_float = refine_peaks(row, ipeaks_int, nwinwidth)[0]

            # self.pintarGrafica(refine_peaks(row, ipeaks_int, nwinwidth)[0] - refinePeaks_spectrum(row, ipeaks_int, nwinwidth))

            fpeaks[idx] = []
            for peak, peak_f in zip(ipeaks_int, ipeaks_float):
                try:
                    sl = numina.array.utils.slice_create(peak, lwidth)
                    rel_peak = peak - sl.start
                    qslit = row[sl]
                    peak_val, fwhm = fmod.compute_fwhm_1d_simple(qslit, rel_peak)
                    peak_on_trace = the_pol(peak)
                    fpeaks[idx].append((peak_f, peak_on_trace, fwhm))
                except ValueError as error:
                    _logger.warning('Error %s computing FWHM in fiber %d', error, idx + 1)
                except IndexError as error:
                    _logger.warning('Error %s computing FWHM in fiber %d', error, idx + 1)
            _logger.debug('found %d peaks in fiber %d', len(fpeaks[idx]), idx)
        return fpeaks

    def pintarGrafica(self, diferencia_final):
        import matplotlib.pyplot as plt

        fig = plt.figure(1)
        ax = fig.add_subplot(111)

        ejeX = numpy.arange(len(diferencia_final))
        ax.plot(ejeX, diferencia_final, label="0")
        # lgd = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=4,
        #                 mode="expand")
        lgd = ax.legend(bbox_to_anchor=(0., 1.02, 1., .102),
                        loc='upper center', ncol=4, mode="expand",
                        borderaxespad=0.)
        handles, labels = ax.get_legend_handles_labels()

        fig.savefig('diferencia.eps', format='eps', dpi=1500,
                    bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.draw()
        plt.show()

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
            _logger.debug('Using %d fiber in reference image', fiberid)
            ref = numpy.array(base[fiberid])
            if ref.size == 0:
                _logger.debug('Reference fiber has no lines, skip')
                continue

            _logger.debug('Positions are %s', ref)
            outref = len(ref)
            savelines = {}
            line_fibers[fiberid] = savelines
            for i in range(outref):
                savelines[i] = {}
                savelines[i]['basedata'] = {}
                savelines[i]['centers'] = []
            _logger.debug('Create kd-tree in fiber %d', fiberid)
            kdtree = cKDTree(ref[:, :2])
            for i in range(ntotal):
                if i == center:
                    for j in range(outref):
                        savelines[j]['basedata']['coordinates'] = tuple(
                            ref[j, :2])
                        savelines[j]['centers'].append(ref[j, 2])
                    continue
                _logger.debug('Matching lines in fiber %d in image # %d', fiberid, i)
                comp = numpy.array(data[i][fiberid])

                if comp.size == 0:
                    _logger.debug('No lines in fiber %d in image # %d', fiberid, i)
                    continue
                else:
                    _logger.debug('Using %d lines in fiber %d in image # %d', comp.size, fiberid, i)

                qdis, qidx = kdtree.query(comp[:, :2],
                                          distance_upper_bound=maxdis)
                for compidx, lidx in enumerate(qidx):
                    if lidx < outref:
                        savelines[lidx]['centers'].append(comp[compidx, 2])

            remove_groups = []

            for ir in savelines:
                if len(savelines[ir]['centers']) != ntotal:
                    remove_groups.append(ir)

            for ir in remove_groups:
                _logger.debug('remove group of lines %d in fiber %d', ir,
                              fiberid)
                del savelines[ir]

        return line_fibers

    def reorder_and_fit(self, line_fibers, focii):
        """Fit all the values of FWHM to a 2nd degree polynomial and return minimum."""

        l = sum(len(value) for key, value in line_fibers.items())

        _logger.debug('there are %d groups of lines to fit', l)
        ally = numpy.zeros((len(focii), l))
        final = numpy.zeros((l, 3))
        _logger.debug('focci are %s', focii)
        l = 0
        for i in line_fibers:
            for j in line_fibers[i]:
                ally[:, l] = line_fibers[i][j]['centers']
                final[l, :2] = line_fibers[i][j]['basedata']['coordinates']
                l += 1

        _logger.debug('line widths are %s', ally)
        res = numpy.polyfit(focii, ally, deg=2)
        _logger.debug('fitting to deg 2 polynomial, done')
        best = -res[1] / (2 * res[0])
        final[:, 2] = best

        return final


def get_valid_traces(tracemap):
    valid = []
    for trace in tracemap:
        fitparms = trace['fitparms']
        if len(fitparms) != 0:
            valid.append(trace['fibid'])
    return valid
