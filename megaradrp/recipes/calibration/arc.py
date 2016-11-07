#
# Copyright 2015-2016 Universidad Complutense de Madrid
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

"""Arc Calibration Recipe for Megara"""

from __future__ import division, print_function

import traceback

import numpy
from astropy.io import fits
from scipy.interpolate import interp1d

from numina.core import Requirement, Product, Parameter, DataFrameType
from numina.core.requirements import ObservationResultRequirement
from numina.core.products import LinesCatalog
from numina.array.wavecalib.arccalibration import arccalibration_direct
from numina.array.wavecalib.arccalibration import fit_list_of_wvfeatures
from numina.array.wavecalib.arccalibration import gen_triplets_master
from numina.array.wavecalib.arccalibration import robust_std
from numina.array.peaks.peakdet import find_peaks_indexes, refine_peaks
from numina.array.peaks.detrend import detrend
from numina.flow import SerialFlow
from numina.array import combine

from skimage.feature import peak_local_max

from megaradrp.processing.combine import basic_processing_with_combination
from megaradrp.processing.aperture import ApertureExtractor
from megaradrp.processing.fiberflat import Splitter, FlipLR
from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.products import WavelengthCalibration
from megaradrp.products.wavecalibration import FiberSolutionArcCalibration
import megaradrp.requirements as reqs

from megaradrp.instrument import vph_thr_arc


class ArcCalibrationRecipe(MegaraBaseRecipe):
    """Process ARC images and create WL_CALIBRATION."""

    # Requirements
    obresult = ObservationResultRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_bpm = reqs.MasterBPMRequirement()
    tracemap = reqs.MasterTraceMapRequirement()
    lines_catalog = Requirement(LinesCatalog, 'Catalog of lines')
    polynomial_degree = Parameter(5, 'Polynomial degree of arc calibration')
    nlines = Parameter(20, "Use the 'nlines' brigthest lines of the catalog")
    # Products
    arc_image = Product(DataFrameType)
    arc_rss = Product(DataFrameType)
    master_wlcalib = Product(WavelengthCalibration)
    fwhm_image = Product(DataFrameType)

    def __init__(self):
        super(ArcCalibrationRecipe, self).__init__("0.1.0")

    def run(self, rinput):

        flow1 = self.init_filters(rinput, rinput.obresult.configuration.values)
        img = basic_processing_with_combination(rinput, flow1, method=combine.median)
        hdr = img[0].header
        self.set_base_headers(hdr)

        splitter1 = Splitter()
        calibrator_aper = ApertureExtractor(rinput.tracemap, self.datamodel)
        flipcor = FlipLR()

        flow2 = SerialFlow([splitter1, calibrator_aper, flipcor])

        reduced_rss = flow2(img)
        reduced2d = splitter1.out

        self.logger.info('extract fibers, %i', len(rinput.tracemap.contents))

        current_vph = rinput.obresult.tags['vph']

        if reduced2d[0].header['INSTRUME'] == 'MEGARA':
            threshold = vph_thr_arc['default'][current_vph]['threshold']
            min_distance = vph_thr_arc['default'][current_vph]['min_distance']
        else:
            raise ValueError('INSTRUME keyword is %s', reduced2d[0].header['INSTRUME'])

        data_wlcalib, fwhm_image = self.calibrate_wl(reduced_rss[0].data,
                                                     rinput.lines_catalog,
                                                     rinput.polynomial_degree,
                                                     rinput.tracemap,
                                                     rinput.nlines,
                                                     threshold=threshold,
                                                     min_distance=min_distance)

        data_wlcalib.tags = rinput.obresult.tags
        # WL calibration goes here
        return self.create_result(arc_image=reduced2d, arc_rss=reduced_rss,
                                  master_wlcalib=data_wlcalib,
                                  fwhm_image=fwhm_image)

    def calc_fwhm_of_line(self, row, peak_int, lwidth=20):
        """
        Compute FWHM of lines in spectra
        """
        import numina.array.fwhm as fmod

        # FIXME: this could wrap around the image
        qslit = row[peak_int - lwidth:peak_int + lwidth]
        return fmod.compute_fwhm_1d_simple(qslit, lwidth)

    def calibrate_wl(self, rss, lines_catalog, poldeg, tracemap, nlines,
                     threshold=0.27,
                     min_distance=30):

        wv_master = lines_catalog[:, 0]
        ntriplets_master, ratios_master_sorted, triplets_master_sorted_list = \
            gen_triplets_master(wv_master)

        nwinwidth = 5
        error_contador = 0
        missing_fib = 0

        data_wlcalib = WavelengthCalibration(instrument='MEGARA')

        for trace in tracemap.contents:
            fibid = trace.fibid
            idx = trace.fibid - 1

            if trace.valid:
                row = rss[idx]

                trend = detrend(row)
                fibdata_detrend = row - trend
                # A fix for LR-V Jun 2016 test images
                # that only have lines there

                row = fibdata_detrend

                self.logger.info('Starting row %d, fibid %d', idx, fibid)
                # find peaks (initial search providing integer numbers)
                ipeaks_int = peak_local_max(row, threshold_rel=threshold,
                                             min_distance=min_distance)[:, 0]

                self.logger.debug('ipeaks_int: %s', ipeaks_int)
                # Filter by flux
                if len(ipeaks_int) <= nlines:
                    # nothing to do
                    ipeaks_int_filtered = ipeaks_int
                else:
                    peak_fluxes = row[ipeaks_int]
                    spositions = peak_fluxes.argsort()
                    # Return last nlines elements
                    ipeaks_int_filtered = ipeaks_int[spositions[-nlines:]]
                    ipeaks_int_filtered.sort()

                self.logger.debug('ipeaks_int_filtered: %s', ipeaks_int_filtered)
                ipeaks_float = refine_peaks(row, ipeaks_int_filtered, nwinwidth)[0]

                # if idx==299:
                if False:
                    import matplotlib.pyplot as plt
                    plt.title('fibid %d' % fibid)
                    plt.plot(row)
                    plt.plot(ipeaks_int, row[ipeaks_int], 'ro', alpha=.9, ms=7,
                             label="ipeaks_int")
                    plt.plot(ipeaks_int_filtered, row[ipeaks_int_filtered], 'ro', alpha=.9, ms=7,
                             label="ipeaks_int_filtered")
                    plt.legend()
                    plt.show()

                # FIXME: xchannel ???
                # This comes from Nico's code, so probably pixels
                # will start in 1
                naxis1 = row.shape[0]
                crpix1 = 1.0
                xchannel = numpy.arange(1, naxis1 + 1)

                finterp_channel = interp1d(range(xchannel.size), xchannel,
                                           kind='linear',
                                           bounds_error=False,
                                           fill_value=0.0)
                xpeaks_refined = finterp_channel(ipeaks_float)

                wv_ini_search = int(
                    lines_catalog[0][0] - 1000)  # initially: 3500
                wv_end_search = int(
                    lines_catalog[-1][0] + 1000)  # initially: 4500

                try:

                    self.logger.info('wv_ini_search %s', wv_ini_search)
                    self.logger.info('wv_end_search %s', wv_end_search)

                    list_of_wvfeatures = arccalibration_direct(
                        wv_master,
                        ntriplets_master,
                        ratios_master_sorted,
                        triplets_master_sorted_list,
                        xpeaks_refined,
                        naxis1,
                        crpix1=crpix1,
                        wv_ini_search=wv_ini_search,
                        wv_end_search=wv_end_search,
                        error_xpos_arc=2.3, # initially: 2.0
                        times_sigma_r=3.0,
                        frac_triplets_for_sum=0.50,
                        times_sigma_theil_sen=10.0,
                        poly_degree_wfit=poldeg,
                        times_sigma_polfilt=10.0,
                        times_sigma_cook=10.0,
                        times_sigma_inclusion=5.0
                    )

                    self.logger.info('Solution for row %d completed', idx)
                    self.logger.info('Fitting solution for row %d', idx)
                    solution_wv = fit_list_of_wvfeatures(
                            list_of_wvfeatures,
                            naxis1_arc=naxis1,
                            crpix1=crpix1,
                            poly_degree_wfit=poldeg,
                            weighted=False,
                            debugplot=0,
                            plot_title=None
                        )

                    self.logger.info('linear crval1, cdelt1: %f %f',
                                     solution_wv.cr_linear.crval,
                                     solution_wv.cr_linear.cdelt)

                    self.logger.info('fitted coefficients %s',
                                     solution_wv.coeff)

                    trace_pol = trace.polynomial
                    # Update feature with measurements of Y coord in original image
                    # Peak and FWHM in RSS
                    for feature in solution_wv.features:
                        # Compute Y
                        feature.ypos = trace_pol(feature.xpos)
                        # FIXME: check here FITS vs PYTHON coordinates, etc
                        peak_int = int(feature.xpos)
                        try:
                            peak, fwhm = self.calc_fwhm_of_line(row, peak_int, lwidth=20)
                        except Exception as error:
                            self.logger.error("%s", error)
                            self.logger.error('error in feature %s', feature)
                            # workaround
                            peak = row[peak_int]
                            fwhm = 0.0
                        # I would call this peak instead...
                        feature.peak = peak
                        feature.fwhm = fwhm

                    # if True:
                    #     plt.title('fibid %d' % fibid)
                    #     plt.plot(row)
                    #     plt.plot(ipeaks_int, row[ipeaks_int],'ro', alpha=.9, ms=7, label="ipeaks_int")
                    #     # # plt.plot(ipeaks_int2, row[ipeaks_int2],'gs', alpha=.5 , ms=10)
                    #     plt.legend()
                    #     plt.show()

                    data_wlcalib.contents[fibid] = FiberSolutionArcCalibration(fibid, solution_wv)

                except (ValueError, TypeError, IndexError) as error:
                    self.logger.error("%s", error)
                    self.logger.error('error in row %d, fibid %d', idx, fibid)
                    traceback.print_exc()
                    data_wlcalib.error_fitting.append(fibid)

                    if False:
                        import matplotlib.pyplot as plt
                        plt.title('fibid %d' % fibid)
                        rrow = row[::-1]
                        rpeaks = 4096-ipeaks_int_filtered[::-1]
                        plt.plot(rrow)
                        plt.plot(rpeaks, rrow[rpeaks], 'ro', alpha=.9, ms=7, label="ipeaks_int_filtered")
                        # # plt.plot(ipeaks_int2, row[ipeaks_int2],'gs', alpha=.5 , ms=10)
                        plt.legend()
                        plt.show()
                    error_contador += 1

            else:
                self.logger.info('skipping row %d, fibid %d, not extracted', idx, fibid)
                missing_fib += 1
                data_wlcalib.missing_fibers.append(fibid)

        self.logger.info('Errors in fitting: %s', error_contador)
        self.logger.info('Missing fibers: %s', missing_fib)

        self.logger.info('Generating fwhm_image...')
        image = self.generate_fwhm_image(data_wlcalib.contents)
        fwhm_image = fits.PrimaryHDU(image)
        fwhm_hdulist = fits.HDUList([fwhm_image])

        self.logger.info('End arc calibration')

        return data_wlcalib, fwhm_hdulist

    def generate_fwhm_image(self, solutions):
        from scipy.spatial import cKDTree

        # Each 10 fibers. Comment this to iterate over all fibers instead.
        ##################################################################
        aux = {}
        for key, value in solutions.items():
            if int(key) % 10 == 0:
                aux[key] = value

        solutions = aux
        ##################################################################

        final = []
        for fibercalib in solutions.values():
            for feature in fibercalib.solution.features:
                final.append([feature.xpos, feature.ypos, feature.fwhm])
        final = numpy.asarray(final)

        voronoi_points = numpy.array(final[:, [0, 1]])

        # Cartesian product of X, Y
        x = numpy.arange(2048 * 2)
        y = numpy.arange(2056 * 2)
        test_points = numpy.transpose(
            [numpy.tile(x, len(y)), numpy.repeat(y, len(x))])

        voronoi_kdtree = cKDTree(voronoi_points)

        test_point_dist, test_point_regions = voronoi_kdtree.query(test_points,
                                                                   k=1)
        final_image = test_point_regions.reshape((4112, 4096)).astype('float64')
        final_image[:, :] = final[final_image[:, :].astype('int64'), 2]
        return (final_image)
