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

        self.logger.info('extract fibers, %i', len(rinput.tracemap.tracelist))
        # List of nonextracted fiberids
        fibids_not_traced = [trace.fibid for trace in rinput.tracemap.tracelist if
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

        fits.writeto('rss.fits', rssdata, clobber=True)

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

        # WL calibration goes here
        return self.create_result(arc_image=reduced, arc_rss=rss,
                                  master_wlcalib=data_wlcalib,
                                  fwhm_image=fwhm_image)

    def run_on_image(self, rssdata, tracemap, threshold_rel, min_distance,
                     limit):
        '''
        Extract spectra, find peaks and compute FWHM.
        '''
        import numina.array.fwhm as fmod
        # Extract the polynomials
        # FIXME: a little hackish
        pols = [numpy.poly1d(t.fitparms) for t in tracemap.tracelist]

        nwinwidth = 5
        lwidth = 20
        fpeaks = {}
        for idx in range(0, len(rssdata)):
            # sampling every 10 fibers...
            row = rssdata[idx, :]
            if numpy.any(row):
                the_pol = pols[idx]
                # find peaks
                trend = detrend(row)
                fibdata_detrend = row - trend
                if limit:
                    row = fibdata_detrend[:limit]
                else:
                    row = fibdata_detrend

                self.logger.info('detect peaks in row %d', idx)
                # find peaks (initial search providing integer numbers)
                ipeaks_int = peak_local_max(row, threshold_rel=threshold_rel,
                                            min_distance=min_distance)[:, 0]

                ipeaks_float = refine_peaks(row, ipeaks_int, nwinwidth)[0]

                fpeaks[idx] = []
                try:
                    for peak, peak_f in zip(ipeaks_int, ipeaks_float):
                        qslit = row[peak - lwidth:peak + lwidth]
                        peak_val, fwhm = fmod.compute_fwhm_1d_simple(qslit,
                                                                     lwidth)
                        peak_on_trace = the_pol(peak)
                        fpeaks[idx].append(
                            (peak_f, peak_on_trace, fwhm, peak_val))
                except:
                    fpeaks[idx] = []
                    # _logger.debug('found %d peaks in fiber %d', len(fpeaks[idx]),idx)
            else:
                fpeaks[idx] = []
        return fpeaks

    def calibrate_wl(self, rss, lines_catalog, poldeg, tracemap,
                     times_sigma=50.0, skiptraces=None, threshold=0.27,
                     min_distance=30):
        #
        # read master table (TBM) and generate auxiliary parameters (valid for
        # all the slits) for the wavelength calibration
        lista_solution = []
        lista_xpeaks_refined = []
        if skiptraces is None:
            skiptraces = []
        wv_master = lines_catalog[:, 0]
        ntriplets_master, ratios_master_sorted, triplets_master_sorted_list = \
            gen_triplets_master(wv_master)
        # FIXME: this depends on the spectral and dispersion axes
        nspec = rss.shape[0]
        coeff_table = numpy.zeros((nspec, poldeg + 1))
        # Loop over rows in RSS
        nwinwidth = 5
        error_contador = 0
        missing_fib = 0

        # # limit = 1900 #LR-Z
        # # limit = 2800 #LR-Z sci
        # # limit = 1900 #LR-V eng ThrNe
        # limit = 0 #LR-V sci

        limit = 0  # LR-R
        dict_of_solution_wv = {}

        # FIXME: hardcoded
        flux_limit = 600000

        for idx, row in enumerate(rss):

            fibid = idx + 1
            if idx==1:
                numpy.savetxt('fiber1.txt', row)
            if fibid not in skiptraces and numpy.any(row):

                trend = detrend(row)
                fibdata_detrend = row - trend
                # A fix for LR-V Jun 2016 test images
                # that only have lines there
                if limit:
                    row = fibdata_detrend[:limit]  # LR-Z
                else:
                    row = fibdata_detrend

                self.logger.info('Starting row %d, fibid %d', idx, fibid)
                # find peaks (initial search providing integer numbers)
                ipeaks_int1 = peak_local_max(row, threshold_rel=threshold,
                                             min_distance=min_distance)[:, 0]
                # filter by flux
                self.logger.info('Filtering peaks over %5.0f', flux_limit)
                ipeaks_vals = row[ipeaks_int1]
                mask = ipeaks_vals < flux_limit
                ipeaks_int = ipeaks_int1[mask]
                self.logger.debug('LEN (ipeaks_int): %s', len(ipeaks_int))
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
                lista_xpeaks_refined.append(xpeaks_refined)
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
                                     solution_wv.crval1_linear,
                                     solution_wv.cdelt1_linear)

                    self.logger.info('fitted coefficients %s',
                                     solution_wv.coeff)
                    coeff_table[idx] = solution_wv.coeff

                    # if True:
                    #     plt.title('fibid %d' % fibid)
                    #     plt.plot(row)
                    #     plt.plot(ipeaks_int, row[ipeaks_int],'ro', alpha=.9, ms=7, label="ipeaks_int")
                    #     # # plt.plot(ipeaks_int2, row[ipeaks_int2],'gs', alpha=.5 , ms=10)
                    #     plt.legend()
                    #     plt.show()
                    dict_of_solution_wv[fibid] = solution_wv

                except (ValueError, TypeError, IndexError) as error:
                    self.logger.error("%s", error)
                    self.logger.error('error in row %d, fibid %d', idx, fibid)
                    traceback.print_exc()
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
                lista_xpeaks_refined.append(numpy.array([]))

            # coeff_table[idx] = numpy_array_with_coeff
        self.logger.info('Errors in fitting: %s', error_contador)
        self.logger.info('Missing fibers: %s', missing_fib)
        lines_rss_fwhm = self.run_on_image(rss, tracemap, threshold,
                                           min_distance, limit)

        data_wlcalib = WavelengthCalibration(instrument='MEGARA')
        data_wlcalib.wvlist = dict_of_solution_wv
        # data_wlcalib = self.generateJSON(coeff_table, dict_of_solution_wv,
        #                                  lista_xpeaks_refined, poldeg,
        #                                  lines_catalog, lines_rss_fwhm)

        self.logger.info('Generating fwhm_image...')
        image = self.generate_image(lines_rss_fwhm)
        fwhm_image = fits.PrimaryHDU(image)
        fwhm_image = fits.HDUList([fwhm_image])

        self.logger.info('End arc calibration')

        return data_wlcalib, fwhm_image

    def generateJSON(self, coeff_table, list_of_solution_wv,
                     lista_xpeaks_refined, poldeg, lines_catalog,
                     lines_rss_fwhm):
        '''
            Final format of the features field is:{
                  "features": [[<x_position>,
                                calculated_wavelength,
                                <original_lambda>,
                                <original_flux>
                                ]]
        '''

        from numpy.polynomial.polynomial import polyval
        self.logger.info('start JSON generation')
        result = []

        for ind, xpeaks in enumerate(lista_xpeaks_refined):
            features = []
            try:
                if numpy.any(xpeaks) and (list_of_solution_wv[ind] is not None) and len(
                        list_of_solution_wv[ind]) == len(xpeaks):
                    res = polyval(xpeaks, coeff_table[ind])
                    # _logger.info('indice: %s', ind)

                    if numpy.any(res):
                        for aux, elem in enumerate(xpeaks):
                            if len(lines_rss_fwhm[ind]) > aux:
                                feature = {'xpos': xpeaks[aux],
                                           'wavelength': res[aux],
                                           'flux': lines_rss_fwhm[ind][aux][3],
                                           'fwhm': lines_rss_fwhm[ind][aux][2],
                                           'ypos': lines_rss_fwhm[ind][aux][1] + 1}

                                features.append(feature)
            except:
                features = []

            function = {
                'method': 'least squares',
                'order': poldeg,
                'coefficients': coeff_table[ind].tolist() if numpy.any(
                    coeff_table[ind]) else []
            }

            record = {}
            record['aperture'] = {'id': ind + 1,
                                  'features': features,
                                  'function': function}
            result.append(record)

        self.logger.info('end JSON generation')

        return result

    def generate_image(self, lines_rss_fwhm):
        from scipy.spatial import cKDTree

        # Each 10 fibers. Comment this to iterate over all fibers instead.
        ##################################################################
        aux = {}
        for key, value in lines_rss_fwhm.items():
            if int(key) % 10 == 0:
                aux[key] = value

        lines_rss_fwhm = aux
        ##################################################################

        l = sum(len(value) for key, value in lines_rss_fwhm.items())

        final = numpy.zeros((l, 3))
        l = 0
        for key, value in lines_rss_fwhm.items():
            for j in range(len(lines_rss_fwhm[key])):
                final[l, :] = lines_rss_fwhm[key][j][:-1]
                l += 1

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
