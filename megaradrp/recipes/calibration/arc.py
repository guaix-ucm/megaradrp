#
# Copyright 2015-2017 Universidad Complutense de Madrid
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
import collections
from datetime import datetime

import numpy
from astropy.io import fits

from numina.core import Requirement, Product, Parameter, DataFrameType
from numina.core.requirements import ObservationResultRequirement
from numina.core.products import LinesCatalog
from numina.array.display.pause_debugplot import pause_debugplot
from numina.array.display.ximplot import ximplot
from numina.array.wavecalib.arccalibration import arccalibration_direct
from numina.array.wavecalib.arccalibration import fit_list_of_wvfeatures
from numina.array.wavecalib.arccalibration import gen_triplets_master
from numina.array.peaks.peakdet import refine_peaks
from numina.array.peaks.detrend import detrend
from numina.flow import SerialFlow
from numina.array import combine
from skimage.feature import peak_local_max

from megaradrp.types import ProcessedFrame, ProcessedRSS
from megaradrp.processing.combine import basic_processing_with_combination
from megaradrp.processing.aperture import ApertureExtractor
from megaradrp.processing.fiberflat import Splitter, FlipLR
from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.products import WavelengthCalibration
from megaradrp.products.wavecalibration import FiberSolutionArcCalibration
import megaradrp.requirements as reqs

from megaradrp.instrument import vph_thr_arc


class ArcCalibrationRecipe(MegaraBaseRecipe):
    """Provides wavelength calibration information from arc images.

    This recipes process a set of arc images obtained in
    *Arc Calibration* mode and returns the information required
    to perform wavelength calibration and resampling in other recipes,
    in the form of a WavelengthCalibration object. The recipe also returns
    a 2D map of the FWHM of the arc lines used for the calibration,
    the result images up to dark correction, and the result of the processing
    up to aperture extraction.


    See Also
    --------
    megaradrp.products.WavelengthCalibration: description of WavelengthCalibration product
    megaradrp.products.LinesCatalog: description of the catalog of lines
    megaradrp.processing.aperture: aperture extraction
    numina.array.wavecalib.arccalibration: wavelength calibration algorithm
    megaradrp.instrument.configs: instrument configuration

    Notes
    -----
    Images provided in `obresult` are trimmed and corrected from overscan,
    bad pixel mask (if `master_bpm` is not None), bias and dark current
    (if `master_dark` is not None).
    Images thus corrected are the stacked using the median.

    The result of the combination is saved as an intermediate result, named
    'reduced_image.fits'. This combined image is also returned in the field
    `reduced_image` of the recipe result.

    The apertures in the 2D image are extracted, using the information in
    `master_traces`. the result of the extraction is saved as an intermediate
    result named 'reduced_rss.fits'. This RSS is also returned in the field
    `reduced_rss` of the recipe result.

    For each fiber in the reduced RSS, the peaks are detected and sorted
    by peak flux. `nlines` is used to select the brightest peaks. If it is a
    list, then the peaks are divided, by their position, in as many groups as
    elements in the list and `nlines[0]` peaks are selected in the first
    group, `nlines[1]` peaks in the second, etc.

    The selected peaks are matched against the catalog of lines in `lines_catalog`.
    The matched lines, the quality of the match and other relevant information is
    stored in the product WavelengthCalibration object. The wavelength of the matched
    features is fitted to a polynomial
    of degree `polynomial_degree`. The coefficients of the polynomial are
    stored in the final `master_wlcalib` object for each fiber.

    """

    # Requirements
    obresult = ObservationResultRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_bpm = reqs.MasterBPMRequirement()
    master_traces = reqs.MasterTraceMapRequirement()
    lines_catalog = Requirement(LinesCatalog, 'Catalog of lines')
    polynomial_degree = Parameter(5, 'Polynomial degree of arc calibration')
    nlines = Parameter(20, "Use the 'nlines' brigthest lines of the spectrum")
    debug_plot = Parameter(0, 'Save intermediate tracing plots')

    # Products
    reduced_image = Product(ProcessedFrame)
    reduced_rss = Product(ProcessedRSS)
    master_wlcalib = Product(WavelengthCalibration)
    fwhm_image = Product(DataFrameType)

    def run(self, rinput):
        """Execute the recipe.

        Parameters
        ----------
        rinput : ArcCalibrationRecipe.RecipeInput

        Returns
        -------
        ArcCalibrationRecipe.RecipeResult

        """

        debugplot = rinput.debug_plot if self.intermediate_results else 0

        flow1 = self.init_filters(rinput, rinput.obresult.configuration)
        img = basic_processing_with_combination(rinput, flow1, method=combine.median)
        hdr = img[0].header
        self.set_base_headers(hdr)

        self.save_intermediate_img(img, 'reduced_image.fits')

        splitter1 = Splitter()
        calibrator_aper = ApertureExtractor(rinput.master_traces, self.datamodel)
        flipcor = FlipLR()

        flow2 = SerialFlow([splitter1, calibrator_aper, flipcor])

        reduced_rss = flow2(img)
        self.save_intermediate_img(reduced_rss, 'reduced_rss.fits')

        reduced2d = splitter1.out

        self.logger.info('extract fibers, %i', len(rinput.master_traces.contents))

        current_vph = rinput.obresult.tags['vph']
        current_insmode = rinput.obresult.tags['insmode']

        if current_insmode in vph_thr_arc and current_vph in vph_thr_arc[current_insmode]:
            threshold = vph_thr_arc[current_insmode][current_vph]['threshold']
            min_distance = vph_thr_arc[current_insmode][current_vph]['min_distance']
            self.logger.info('rel threshold for %s is %4.2f', current_vph, threshold)
        else:
            threshold = 0.02
            min_distance = 10.0
            self.logger.info('rel threshold not defined for %s, using %4.2f', current_vph, threshold)
            self.logger.info('min dist not defined for %s, using %4.2f', current_vph, min_distance)

        if isinstance(rinput.nlines, collections.Iterable):
            nlines = rinput.nlines
        else:
            nlines = [rinput.nlines]

        data_wlcalib, fwhm_image = self.calibrate_wl(reduced_rss[0].data,
                                                     rinput.lines_catalog,
                                                     rinput.polynomial_degree,
                                                     rinput.master_traces,
                                                     nlines,
                                                     threshold=threshold,
                                                     min_distance=min_distance,
                                                     debugplot=debugplot)

        data_wlcalib.tags = rinput.obresult.tags
        final = data_wlcalib
        final.meta_info['creation_date'] = datetime.utcnow().isoformat()
        final.meta_info['mode_name'] = self.mode
        final.meta_info['instrument_name'] = self.instrument
        final.meta_info['recipe_name'] = self.__class__.__name__
        final.meta_info['recipe_version'] = self.__version__
        final.meta_info['origin'] = {}
        final.meta_info['origin']['block_uuid'] = reduced2d[0].header.get('BLCKUUID', "UNKNOWN")
        final.meta_info['origin']['insconf_uuid'] = reduced2d[0].header.get('INSCONF', "UNKNOWN")
        final.meta_info['origin']['date_obs'] = reduced2d[0].header['DATE-OBS']

        # FIXME: redundant
        cdata = []
        for frame in rinput.obresult.frames:
            hdulist = frame.open()
            fname = self.datamodel.get_imgid(hdulist)
            cdata.append(fname)

        final.meta_info['origin']['frames'] = cdata
        # WL calibration goes here
        return self.create_result(reduced_image=reduced2d, reduced_rss=reduced_rss,
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
                     min_distance=30,
                     debugplot=0):

        wv_master = lines_catalog[:, 0]
        ntriplets_master, ratios_master_sorted, triplets_master_sorted_list = \
            gen_triplets_master(wv_master)

        nwinwidth = 5
        error_contador = 0
        missing_fib = 0

        data_wlcalib = WavelengthCalibration(instrument='MEGARA')
        data_wlcalib.total_fibers = tracemap.total_fibers
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

                self.logger.debug('ipeaks_int.........: %s', ipeaks_int)

                # Filter by flux, selecting a maximum number of brightest
                # lines in each region
                region_size = (len(row)-1)/(len(nlines))
                ipeaks_int_filtered = numpy.array([], dtype=int)
                for iregion, nlines_in_region in enumerate(nlines):
                    if nlines_in_region > 0:
                        imin = int(iregion * region_size)
                        imax = int((iregion + 1) * region_size)
                        if iregion > 0:
                            imin += 1
                        ipeaks_int_region = ipeaks_int[numpy.logical_and(
                            ipeaks_int >= imin, ipeaks_int <= imax
                        )]
                        if len(ipeaks_int_region) > 0:
                            peak_fluxes = row[ipeaks_int_region]
                            spos = peak_fluxes.argsort()
                            ipeaks_tmp = ipeaks_int_region[
                                spos[-nlines_in_region:]
                            ]
                            ipeaks_tmp.sort()  # in-place sort
                            self.logger.debug('ipeaks_in_region...: %s',
                                              ipeaks_tmp)
                            ipeaks_int_filtered=numpy.concatenate(
                                (ipeaks_int_filtered, ipeaks_tmp)
                            )

                self.logger.debug('ipeaks_int_filtered: %s', ipeaks_int_filtered)
                ipeaks_float = refine_peaks(row, ipeaks_int_filtered, nwinwidth)[0]

                # if 43 <= fibid <= 44:
                #    debugplot = 12
                #else:
                #    debugplot = 0
                if abs(debugplot) > 10:
                    print('ipeaks_float:\n', ipeaks_float)
                    ax = ximplot(row, show=False, debugplot=debugplot)
                    ax.set_title('fibid %d' % fibid)
                    ax.plot(row)
                    ax.plot(ipeaks_int, row[ipeaks_int], 'b+', alpha=.9, ms=7,
                            label="ipeaks_int")
                    ax.plot(ipeaks_int_filtered, row[ipeaks_int_filtered],
                            'ro', alpha=.9, ms=7, label="ipeaks_int_filtered")
                    ax.plot(ipeaks_float, row[ipeaks_int_filtered], 'go',
                            alpha=.9, ms=7, label="ipeaks_float")
                    ax.legend()
                    pause_debugplot(debugplot, pltshow=True)

                # FIXME: xchannel ???
                naxis1 = row.shape[0]
                crpix1 = 1.0
                # This comes from Nico's code, so probably pixels
                # will start in 1
                # xchannel = numpy.arange(1, naxis1 + 1)
                #
                # finterp_channel = interp1d(range(xchannel.size), xchannel,
                #                            kind='linear',
                #                            bounds_error=False,
                #                            fill_value=0.0)
                # xpeaks_refined = finterp_channel(ipeaks_float)
                xpeaks_refined = ipeaks_float + 1.0

                wv_range_catalog = lines_catalog[-1][0] - lines_catalog[0][0]
                delta_wv = 0.20 * wv_range_catalog
                wv_ini_search = int(lines_catalog[0][0] - delta_wv)
                wv_end_search = int(lines_catalog[-1][0] + delta_wv)

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
                        error_xpos_arc=3.0, # initially: 2.0
                        times_sigma_r=3.0,
                        frac_triplets_for_sum=0.50,
                        times_sigma_theil_sen=10.0,
                        poly_degree_wfit=poldeg,
                        times_sigma_polfilt=10.0,
                        times_sigma_cook=10.0,
                        times_sigma_inclusion=5.0,
                        debugplot=debugplot
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

                    new = FiberSolutionArcCalibration(fibid, solution_wv)
                    data_wlcalib.contents.append(new)

                except (ValueError, TypeError, IndexError) as error:
                    self.logger.error("%s", error)
                    self.logger.error('error in row %d, fibid %d', idx, fibid)
                    traceback.print_exc()
                    data_wlcalib.error_fitting.append(fibid)

                    if abs(debugplot) > 10:
                        rrow = row[::-1]
                        rpeaks = 4096-ipeaks_int_filtered[::-1]
                        ax = ximplot(rrow, show=False, debugplot=debugplot)
                        ax.set_title('fibid %d' % fibid)
                        ax.plot(rpeaks, rrow[rpeaks], 'ro', alpha=.9, ms=7,
                                label="ipeaks_int_filtered")
                        # ax.plot(ipeaks_int2, row[ipeaks_int2],'gs',
                        #         alpha=.5 , ms=10)
                        ax.legend()
                        pause_debugplot(debugplot, pltshow=True)
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
        aux = []
        for value in solutions:
            if value.fibid % 10 == 0:
                aux.append(value)

        ##################################################################

        final = []
        for fibercalib in aux:
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
        final_image = test_point_regions.reshape((4112, 4096)).astype('float32')
        final_image[:, :] = final[final_image[:, :].astype('int16'), 2]
        return (final_image)
