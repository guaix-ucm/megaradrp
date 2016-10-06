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
from skimage.feature import peak_local_max

from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.products import WavelengthCalibration
from megaradrp.products.wavecalibration import FiberSolutionArcCalibration
import megaradrp.requirements as reqs
from megaradrp.core.processing import apextract_tracemap_2


# FIXME: hardcoded numbers
vph_thr = {'default':{'LR-I':{'min_distance':10,
                              'threshold':0.06},
                      'LR-R':{'min_distance':10,
                              'threshold':0.20},
                      'LR-V': {'min_distance':30,
                               'threshold':0.19},
                      'LR-Z': {'min_distance':60,
                               'threshold':0.02},
                      'LR-U':{'min_distance':10,
                              'threshold': 0.02,}
                      },
}


class ArcCalibrationRecipe(MegaraBaseRecipe):
    """Process ARC images and create WL_CALIBRATION."""

    # Requirements
    obresult = ObservationResultRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    # master_bpm = reqs.MasterBPMRequirement()
    tracemap = reqs.MasterTraceMapRequirement()
    lines_catalog = Requirement(LinesCatalog, 'Catalog of lines')
    polynomial_degree = Parameter(5, 'Polynomial degree of arc calibration')
    # Products
    arc_image = Product(DataFrameType)
    arc_rss = Product(DataFrameType)
    master_wlcalib = Product(WavelengthCalibration)
    fwhm_image = Product(DataFrameType)

    def __init__(self):
        super(ArcCalibrationRecipe, self).__init__("0.1.0")

    def run(self, rinput):
        # Basic processing
        parameters = self.get_parameters(rinput)
        reduced = self.bias_process_common(rinput.obresult, parameters)

        self.logger.info('extract fibers, %i', len(rinput.tracemap.contents))
        # List of nonextracted fiberids
        fibids_not_traced = [trace.fibid for trace in rinput.tracemap.contents if
                             not trace.fitparms]
        self.logger.info('not traced fibers, %i', len(fibids_not_traced))

        # rssdata = apextract_tracemap(reduced[0].data, rinput.tracemap)
        rssdata = apextract_tracemap_2(reduced[0].data, rinput.tracemap)
        rssdata = numpy.fliplr(rssdata)

        rsshdu = fits.PrimaryHDU(rssdata, header=reduced[0].header)
        header_list = self.getHeaderList(
            [reduced, rinput.obresult.images[0].open()])
        rss = fits.HDUList([rsshdu] + header_list)

        self.logger.info('extracted %i fibers', rssdata.shape[0])

        # Skip any other inputs for the moment
        # data_wlcalib, fwhm_image = self.calibrate_wl(rssdata, rinput.lines_catalog,
        #                                  rinput.polynomial_degree, rinput.tracemap)

        current_vph = rinput.obresult.tags['vph']

        if reduced[0].header['INSTRUME'] == 'MEGARA':
            threshold = vph_thr['default'][current_vph]['threshold']
            min_distance = vph_thr['default'][current_vph]['min_distance']
        else:
            raise ValueError('INSTRUME keyword is %s', reduced[0].header['INSTRUME'])

        data_wlcalib, fwhm_image = self.calibrate_wl(rssdata,
                                                     rinput.lines_catalog,
                                                     rinput.polynomial_degree,
                                                     rinput.tracemap,
                                                     skiptraces=fibids_not_traced,
                                                     threshold=threshold,
                                                     min_distance=min_distance)

        data_wlcalib.tags = rinput.obresult.tags
        # WL calibration goes here
        return self.create_result(arc_image=reduced, arc_rss=rss,
                                  master_wlcalib=data_wlcalib,
                                  fwhm_image=fwhm_image)

    def calc_fwhm_of_line(self, row, peak_int, lwidth=20):
        '''
        Compute FWHM of lines in spectra
        '''
        import numina.array.fwhm as fmod

        # FIXME: this could wrap around the image
        qslit = row[peak_int - lwidth:peak_int + lwidth]
        return fmod.compute_fwhm_1d_simple(qslit, lwidth)

    def calibrate_wl(self, rss, lines_catalog, poldeg, tracemap,
                     times_sigma=50.0, skiptraces=None, threshold=0.27,
                     min_distance=30):

        if skiptraces is None:
            skiptraces = []
        wv_master = lines_catalog[:, 0]
        ntriplets_master, ratios_master_sorted, triplets_master_sorted_list = \
            gen_triplets_master(wv_master)

        nwinwidth = 5
        error_contador = 0
        missing_fib = 0

        # FIXME: make trace map use new polynomials instead of poly1d
        trace_pols = [numpy.poly1d(t.fitparms) for t in tracemap.contents]

        data_wlcalib = WavelengthCalibration(instrument='MEGARA')

        for idx, row in enumerate(rss):

            fibid = idx + 1

            if fibid not in skiptraces and numpy.any(row):

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
                ipeaks_float = refine_peaks(row, ipeaks_int, nwinwidth)[0]

                # if idx==299:
                if False:
                    import matplotlib.pyplot as plt
                    plt.title('fibid %d' % fibid)
                    plt.plot(row)
                    plt.plot(ipeaks_int, row[ipeaks_int], 'ro', alpha=.9, ms=7,
                             label="ipeaks_int")
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


                    trace_pol = trace_pols[idx]
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
                        rpeaks = 4096-ipeaks_int[::-1]
                        plt.plot(rrow)
                        plt.plot(rpeaks, rrow[rpeaks], 'ro', alpha=.9, ms=7, label="ipeaks_int")
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
