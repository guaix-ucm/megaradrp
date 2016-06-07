#
# Copyright 2015 Universidad Complutense de Madrid
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

import logging

import numpy
from astropy.io import fits
from scipy.interpolate import interp1d
from numina.core import Requirement, Product, Parameter, DataFrameType
from numina.core.requirements import ObservationResultRequirement
from numina.core.products import LinesCatalog
from numina.array.wavecal.arccalibration import arccalibration_direct
from numina.array.wavecal.arccalibration import fit_solution
from numina.array.wavecal.arccalibration import gen_triplets_master
from numina.array.wavecal.statsummary import sigmaG
from numina.array.peaks.peakdet import find_peaks_indexes, refine_peaks

from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.products import TraceMap, WavelengthCalibration, JSONstorage
from megaradrp.requirements import MasterBiasRequirement, MasterBPMRequirement
from megaradrp.requirements import MasterDarkRequirement
from megaradrp.core.processing import apextract_tracemap
from megaradrp.core.processing import apextract_tracemap_2
from skimage.feature import peak_local_max

_logger = logging.getLogger('numina.recipes.megara')


class ArcCalibrationRecipe(MegaraBaseRecipe):
    """Process ARC images and create WL_CALIBRATION."""

    # Requirements
    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_bpm = MasterBPMRequirement()
    tracemap = Requirement(TraceMap, 'Trace information of the Apertures')
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

        _logger.info('extract fibers')
        _logger.info('extract fibers, %i', len(rinput.tracemap))
        # rssdata = apextract_tracemap(reduced[0].data, rinput.tracemap)
        rssdata = apextract_tracemap_2(reduced[0].data, rinput.tracemap)
        rsshdu = fits.PrimaryHDU(rssdata, header=reduced[0].header)
        header_list = self.getHeaderList(
            [reduced, rinput.obresult.images[0].open()])
        rss = fits.HDUList([rsshdu] + header_list)

        _logger.info('extracted %i fibers', rssdata.shape[0])

        fits.writeto('rss.fits', rssdata, clobber=True)

        # Skip any other inputs for the moment
        # data_wlcalib, fwhm_image = self.calibrate_wl(rssdata, rinput.lines_catalog,
        #                                  rinput.polynomial_degree, rinput.tracemap)
        
        data_wlcalib, fwhm_image = self.calibrate_wl2(rssdata, rinput.lines_catalog,
                                        rinput.polynomial_degree, rinput.tracemap)
        
        # WL calibration goes here
        return self.create_result(arc_image=reduced, arc_rss=rss,
                                  master_wlcalib=data_wlcalib, fwhm_image=fwhm_image)

    def run_on_image(self, rssdata, tracemap):
            '''
            Extract spectra, find peaks and compute FWHM.
            '''
            import numina.array.fwhm as fmod
            # Extract the polynomials
            # FIXME: a little hackish
            pols = [numpy.poly1d(t['fitparms']) for t in tracemap]

            nwinwidth = 5
            lwidth = 20
            fpeaks = {}
            for idx in range(0, len(rssdata)):
                # sampling every 10 fibers...
                row = rssdata[idx, :]

                if numpy.any(row):
                    the_pol = pols[idx]
                    # find peaks
                    trend = self.detrend(row)
                    fibdata_detrend = row - trend
                    row = fibdata_detrend[:1850]

                    # _logger.info('Starting row %d', idx)
                    # find peaks (initial search providing integer numbers)
                    ipeaks_int = peak_local_max(row, threshold_rel=0.02, min_distance=40)[:, 0]

                    ipeaks_float = refine_peaks(row, ipeaks_int, nwinwidth)[0]

                    fpeaks[idx] = []
                    for peak, peak_f in zip(ipeaks_int, ipeaks_float):
                        qslit = row[peak - lwidth:peak + lwidth]
                        peak_val, fwhm = fmod.compute_fwhm_1d_simple(qslit, lwidth)
                        peak_on_trace = the_pol(peak)
                        fpeaks[idx].append((peak_f, peak_on_trace, fwhm, peak_val))
                    # _logger.debug('found %d peaks in fiber %d', len(fpeaks[idx]),idx)
                else:
                    fpeaks[idx] = []
            return fpeaks


    def detrend(self, m, deg=5, tol=1e-3):
        import numpy as np
        nloop = 10
        xx = np.arange(len(m))
        ss = m.copy()
        m2 = ss
        pol = np.ones((deg+1,))
        for _ in range(nloop):
            pol_new = np.polyfit(xx, ss, deg)
            pol2 = numpy.linalg.norm(pol)
            pol1 = numpy.linalg.norm(pol - pol_new)
            if pol1 / pol2 < tol:
                break
            pol = pol_new
            m2 = np.polyval(pol, xx)
            ss = np.minimum(ss, m2)
        return m2

    def calibrate_wl2(self, rss, lines_catalog, poldeg, tracemap, times_sigma=50.0):
        #
        # read master table (TBM) and generate auxiliary parameters (valid for
        # all the slits) for the wavelength calibration
        lista_solution = []
        lista_xpeaks_refined = []

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
        for idx, row in enumerate(rss):
            solution = []
            if numpy.any(row):
                try:
                    import matplotlib.pyplot as plt
                    from skimage.filters import threshold_otsu

                    # row = row[:2250] #funciona con idx=299 threshold_rel=0.13, min_distance=50
                    # row = row[:1860]
                    # aux = numpy.ones((2850))
                    # aux[:1350] = row[:1350]
                    # aux[1350:] = row[1450:2950]
                    # row = aux

                    trend = self.detrend(row)
                    fibdata_detrend = row - trend
                    row = fibdata_detrend[:1850]

                    _logger.info('Starting row %d', idx)
                    # find peaks (initial search providing integer numbers)
                    ipeaks_int = peak_local_max(row, threshold_rel=0.02, min_distance=40)[:, 0]
                    logging.warning('ipeaks_int: %s', ipeaks_int)
                    ipeaks_float = refine_peaks(row, ipeaks_int, nwinwidth)[0]

                    # plt.plot(row)
                    # plt.plot(ipeaks_int, row[ipeaks_int],'ro', alpha=.9, ms=7, label="ipeaks_int")
                    # # plt.plot(ipeaks_int2, row[ipeaks_int2],'gs', alpha=.5 , ms=10)
                    # plt.legend()
                    # plt.show()

                    # FIXME: xchannel ???
                    # This comes from Nico's code, so probably pixels
                    # will start in 1
                    naxis1 = row.shape[0]
                    xchannel = numpy.arange(1, naxis1 + 1)

                    finterp_channel = interp1d(range(xchannel.size), xchannel,
                                               kind='linear')
                    xpeaks_refined = finterp_channel(ipeaks_float)
                    lista_xpeaks_refined.append(xpeaks_refined)
                    wv_ini_search = int(lines_catalog[0][0]-1000) # initially: 3500
                    wv_end_search = int(lines_catalog[-1][0]+1000) #initially: 4500

                    _logger.info('wv_ini_search %s', wv_ini_search)
                    _logger.info('wv_end_search %s', wv_end_search)

                    try:
                        solution = arccalibration_direct(wv_master,
                                                         ntriplets_master,
                                                         ratios_master_sorted,
                                                         triplets_master_sorted_list,
                                                         xpeaks_refined,
                                                         naxis1,
                                                         wv_ini_search=wv_ini_search,
                                                         wv_end_search=wv_end_search,
                                                         error_xpos_arc=0.3, #initially: 2.0
                                                         times_sigma_r=3.0,
                                                         frac_triplets_for_sum=0.50,
                                                         times_sigma_theil_sen=10.0,
                                                         poly_degree_wfit=poldeg,
                                                         times_sigma_polfilt=10.0,
                                                         times_sigma_inclusion=5.0)

                        _logger.info('Solution for row %d completed', idx)
                        _logger.info('Fitting solution for row %d', idx)
                        numpy_array_with_coeff, crval1_approx, cdelt1_approx = \
                            fit_solution(wv_master,
                                         xpeaks_refined,
                                         solution,
                                         poly_degree_wfit=poldeg,
                                         weighted=False)

                        _logger.info('approximate crval1, cdelt1: %f %f',
                                     crval1_approx,
                                     cdelt1_approx)

                        _logger.info('fitted coefficients %s', numpy_array_with_coeff)
                        coeff_table[idx] = numpy_array_with_coeff
                    except TypeError as error:
                        _logger.error("%s", error)
                except ValueError as error:
                    _logger.error("%s", error)
                    _logger.error('Erro en Fibra: %s', idx)
                    error_contador += 1
            else:
                missing_fib += 1
                lista_xpeaks_refined.append(numpy.array([]))

            _logger.error('ERRORES: %s', error_contador)
            _logger.error('Perdidas: %s', missing_fib)

            lista_solution.append(solution)
            # coeff_table[idx] = numpy_array_with_coeff


        lines_rss_fwhm = self.run_on_image(rss, tracemap)
        data_wlcalib = self.generateJSON(coeff_table, lista_solution,
                                         lista_xpeaks_refined, poldeg,
                                         lines_catalog, lines_rss_fwhm)

        _logger.info('Generating fwhm_image...')
        image = self.generate_image(lines_rss_fwhm)
        fwhm_image = fits.PrimaryHDU(image)
        fwhm_image = fits.HDUList([fwhm_image])

        _logger.info('Fin')

        return data_wlcalib, fwhm_image


    def calibrate_wl(self, rss, lines_catalog, poldeg, tracemap, times_sigma=50.0):
        #
        # read master table (TBM) and generate auxiliary parameters (valid for
        # all the slits) for the wavelength calibration
        lista_solution = []
        lista_xpeaks_refined = []

        wv_master = lines_catalog[:, 0]
        ntriplets_master, ratios_master_sorted, triplets_master_sorted_list = \
            gen_triplets_master(wv_master)
        # FIXME: this depends on the spectral and dispersion axes
        nspec = rss.shape[0]
        coeff_table = numpy.zeros((nspec, poldeg + 1))
        # Loop over rows in RSS
        nwinwidth = 5
        for idx, row in enumerate(rss):
            _logger.info('Starting row %d', idx)
            # find peaks (initial search providing integer numbers)
            threshold = numpy.median(row) + times_sigma * sigmaG(row)

            ipeaks_int = find_peaks_indexes(row, nwinwidth, threshold)
            ipeaks_float = refine_peaks(row, ipeaks_int, nwinwidth)[0]

            # define interpolation function and interpolate the refined peak
            # location, passing from index number (within the row array)
            # to channel number (note that this step takes care of the fact
            # that the extracted spectrum may correspond to a subregion in the
            # spectral direction)

            # FIXME: xchannel ???
            # This comes from Nico's code, so probably pixels
            # will start in 1
            naxis1 = row.shape[0]
            xchannel = numpy.arange(1, naxis1 + 1)

            finterp_channel = interp1d(range(xchannel.size), xchannel,
                                       kind='linear')
            xpeaks_refined = finterp_channel(ipeaks_float)
            lista_xpeaks_refined.append(xpeaks_refined)
            wv_ini_search = int(lines_catalog[0][0] - 1000)  # initially: 3500
            wv_end_search = int(lines_catalog[-1][0] + 1000)  # initially: 4500

            _logger.info('wv_ini_search %s', wv_ini_search)
            _logger.info('wv_end_search %s', wv_end_search)

            try:
                solution = arccalibration_direct(wv_master,
                                                 ntriplets_master,
                                                 ratios_master_sorted,
                                                 triplets_master_sorted_list,
                                                 xpeaks_refined,
                                                 naxis1,
                                                 wv_ini_search=wv_ini_search,
                                                 wv_end_search=wv_end_search,
                                                 error_xpos_arc=0.3, #initially: 2.0
                                                 times_sigma_r=3.0,
                                                 frac_triplets_for_sum=0.50,
                                                 times_sigma_theil_sen=10.0,
                                                 poly_degree_wfit=poldeg,
                                                 times_sigma_polfilt=10.0,
                                                 times_sigma_inclusion=5.0)

                _logger.info('Solution for row %d completed', idx)
                _logger.info('Fitting solution for row %d', idx)
                numpy_array_with_coeff, crval1_approx, cdelt1_approx = \
                    fit_solution(wv_master,
                                 xpeaks_refined,
                                 solution,
                                 poly_degree_wfit=poldeg,
                                 weighted=False)

                _logger.info('approximate crval1, cdelt1: %f %f',
                             crval1_approx,
                             cdelt1_approx)

                _logger.info('fitted coefficients %s', numpy_array_with_coeff)
                lista_solution.append(solution)
                coeff_table[idx] = numpy_array_with_coeff

            except TypeError as error:
                _logger.error("%s", error)

        lines_rss_fwhm = self.run_on_image(rss, tracemap)
        data_wlcalib = self.generateJSON(coeff_table, lista_solution,
                                         lista_xpeaks_refined, poldeg,
                                         lines_catalog, lines_rss_fwhm)

        _logger.info('Generating fwhm_image...')
        image = self.generate_image(lines_rss_fwhm)
        fwhm_image = fits.PrimaryHDU(image)
        fwhm_image = fits.HDUList([fwhm_image])

        _logger.info('Fin')

        return data_wlcalib, fwhm_image

    def generateJSON(self, coeff_table, lista_solution,
                     lista_xpeaks_refined, poldeg, lines_catalog, lines_rss_fwhm):
        '''
            Final format of the features field is:{
                  "features": [[<x_position>,
                                calculated_wavelength,
                                <original_lambda>,
                                <original_flux>
                                ]]
        '''

        from numpy.polynomial.polynomial import polyval
        _logger.info('start JSON generation')
        result = []
        
        for ind, xpeaks in enumerate(lista_xpeaks_refined):
            features = []
            
            if numpy.any(xpeaks) and lista_solution[ind]:
                res = polyval(xpeaks, coeff_table[ind])
                # _logger.info('indice: %s', ind)
                
                if numpy.any(res):
                    for aux, elem in enumerate(xpeaks):
                        feature = {'xpos':None,
                                   'wavelength':None,
                                   'reference':None,
                                   'flux':None,
                                   'category':None
                          }

                        feature['xpos'] = xpeaks[aux]
                        feature['wavelength'] = res[aux]
                        feature['reference'] = lines_catalog[aux][0]
                        feature['flux'] = lines_rss_fwhm[ind][aux][3]
                        feature['category'] = lista_solution[ind][aux]['type']
                        feature['fwhm'] = lines_rss_fwhm[ind][aux][2]
                        feature['ypos'] = lines_rss_fwhm[ind][aux][1] + 1## Fits coordinate.
                        features.append(feature)

            function = {
                'method':'least squares',
                'order':poldeg,
                'coefficients': coeff_table[ind].tolist() if numpy.any(coeff_table[ind]) else []
            }

            record = {}
            record['aperture'] = {'id': ind + 1,
                                  'features': features,
                                  'function': function}
            result.append(record)

        _logger.info('end JSON generation')

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

        test_points = numpy.transpose([numpy.tile(x, len(y)), numpy.repeat(y, len(x))])

        voronoi_kdtree = cKDTree(voronoi_points)

        test_point_dist, test_point_regions = voronoi_kdtree.query(test_points,k=1)
        final_image = test_point_regions.reshape((4112, 4096)).astype('float64')
        final_image[:, :] = final[final_image[:, :].astype('int64'), 2]
        return (final_image)
