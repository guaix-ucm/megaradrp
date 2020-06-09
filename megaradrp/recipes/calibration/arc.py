#
# Copyright 2015-2019 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Arc Calibration Recipe for Megara"""

from __future__ import division, print_function

import traceback
from datetime import datetime

import numpy
from astropy.io import fits
from copy import deepcopy
import errno
import os

from numina.core import Requirement, Result, Parameter, DataFrameType
from numina.core.requirements import ObservationResultRequirement
from numina.array.display.polfit_residuals import polfit_residuals
from numina.array.display.polfit_residuals import \
    polfit_residuals_with_sigma_rejection
from numina.array.display.ximplotxy import ximplotxy
from numina.array.wavecalib.__main__ import find_fxpeaks
from numina.array.wavecalib.arccalibration import arccalibration_direct
from numina.array.wavecalib.arccalibration import fit_list_of_wvfeatures
from numina.array.wavecalib.arccalibration import gen_triplets_master
from numina.array.wavecalib.arccalibration import refine_arccalibration
from numina.array.wavecalib.solutionarc import CrLinear
from numina.array.wavecalib.solutionarc import SolutionArcCalibration
from numina.core.validator import range_validator
from numina.util.flow import SerialFlow
from numina.array import combine

from megaradrp.ntypes import ProcessedFrame, ProcessedRSS
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
    master_apertures = reqs.MasterAperturesRequirement()
    extraction_offset = Parameter([0.0], 'Offset traces for extraction',
                                  accept_scalar=True)
    lines_catalog = reqs.LinesCatalogRequirement()
    polynomial_degree = Parameter(5, 'Polynomial degree of arc calibration',
                                  as_list=True, nelem='+',
                                  validator=range_validator(minval=1)
                                  )
    nlines = Parameter(20, "Use the 'nlines' brigthest lines of the spectrum",
                       as_list=True, nelem='+',
                       validator=range_validator(minval=0))
    debug_plot = Parameter(0, 'Save intermediate tracing plots')
    store_pdf_with_refined_fits = Parameter(
        0,
        description='Store PDF plot with refined fits for each fiber',
    )

    # Results
    reduced_image = Result(ProcessedFrame)
    reduced_rss = Result(ProcessedRSS)
    master_wlcalib = Result(WavelengthCalibration)
    fwhm_image = Result(DataFrameType)

    def run(self, rinput):
        """Execute the recipe.

        Parameters
        ----------
        rinput : ArcCalibrationRecipe.RecipeInput

        Returns
        -------
        ArcCalibrationRecipe.RecipeResult

        """

        self.logger.info('starting arc calibration recipe')

        debugplot = rinput.debug_plot if self.intermediate_results else 0

        obresult = rinput.obresult
        obresult_meta = obresult.metadata_with(self.datamodel)

        flow1 = self.init_filters(rinput, rinput.obresult.configuration)
        img = basic_processing_with_combination(rinput, flow1, method=combine.median)
        hdr = img[0].header
        self.set_base_headers(hdr)

        self.save_intermediate_img(img, 'reduced_image.fits')

        splitter1 = Splitter()
        calibrator_aper = ApertureExtractor(
            rinput.master_apertures,
            self.datamodel,
            offset=rinput.extraction_offset
        )
        flipcor = FlipLR()

        flow2 = SerialFlow([splitter1, calibrator_aper, flipcor])

        reduced_rss = flow2(img)
        self.save_intermediate_img(reduced_rss, 'reduced_rss.fits')

        reduced2d = splitter1.out

        self.logger.info('extract fibers, %i', len(rinput.master_apertures.contents))

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

        # WL calibration goes here
        initial_data_wlcalib, data_wlcalib, fwhm_image = self.calibrate_wl(
            reduced_rss[0].data,
            rinput.lines_catalog,
            rinput.polynomial_degree,
            rinput.master_apertures, rinput.nlines,
            threshold=threshold,
            min_distance=min_distance,
            debugplot=debugplot,
            store_pdf_with_refined_fits=rinput.store_pdf_with_refined_fits
        )

        initial_data_wlcalib.tags = rinput.obresult.tags
        final = initial_data_wlcalib
        final.update_metadata(self)
        final.update_metadata_origin(obresult_meta)

        self.save_structured_as_json(
            initial_data_wlcalib,
            'initial_master_wlcalib.json'
        )

        self.logger.info('end arc calibration recipe')

        if data_wlcalib is None:
            return self.create_result(
                reduced_image=reduced2d,
                reduced_rss=reduced_rss,
                master_wlcalib=initial_data_wlcalib,
                fwhm_image=fwhm_image
            )
        else:
            # copy metadata from initial_master_wlcalib to master_wlcalib
            data_wlcalib.tags = deepcopy(initial_data_wlcalib.tags)
            data_wlcalib.meta_info = deepcopy(initial_data_wlcalib.meta_info)
            return self.create_result(
                reduced_image=reduced2d,
                reduced_rss=reduced_rss,
                master_wlcalib=data_wlcalib,
                fwhm_image=fwhm_image
            )

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
                     debugplot=0,
                     store_pdf_with_refined_fits=0):

        if len(poldeg) == 1:
            poldeg_initial = poldeg[0]
            poldeg_refined = poldeg[0]
        elif len(poldeg) == 2:
            poldeg_initial = poldeg[0]
            poldeg_refined = poldeg[1]
        else:
            raise ValueError("Unexpected list length for polynomial_degree")

        if poldeg_initial > poldeg_refined:
            raise ValueError("Unexpected poldeg_initial (" +
                             str(poldeg_initial) + ") > poldeg_refined (" +
                             str(poldeg_refined) + ")")

        wv_master_all = lines_catalog[:, 0]
        if lines_catalog.shape[1] == 2:  # assume old format
            wv_master = numpy.copy(wv_master_all)
            refine_wv_calibration = False
            if abs(debugplot) >= 10:
                print('wv_master:\n', wv_master)
        elif lines_catalog.shape[1] == 3:  # assume new format
            wv_flag = lines_catalog[:, 1]
            wv_master = wv_master_all[numpy.where(wv_flag == 1)]
            refine_wv_calibration = True
            if abs(debugplot) >= 10:
                print('wv_master:\n', wv_master)
        else:
            raise ValueError('lines_catalog file does not have the expected '
                             'number of columns')

        # minimum and maximum wavelength
        wv_range_catalog = wv_master_all[-1] - wv_master_all[0]
        delta_wv = 0.20 * wv_range_catalog
        wv_ini_search = int(wv_master_all[0] - delta_wv)
        wv_end_search = int(wv_master_all[-1] + delta_wv)
        self.logger.info('wv_ini_search %s', wv_ini_search)
        self.logger.info('wv_end_search %s', wv_end_search)

        ntriplets_master, ratios_master_sorted, triplets_master_sorted_list = \
            gen_triplets_master(wv_master)

        error_contador = 0
        missing_fib = 0
        crpix1 = 1.0
        naxis1 = rss.shape[1]

        plot_tracenumber = []
        plot_npeaksfound = []
        plot_crval1 = []
        plot_cdelt1 = []
        plot_coeff = []

        initial_data_wlcalib = WavelengthCalibration(instrument='MEGARA')
        initial_data_wlcalib.total_fibers = tracemap.total_fibers
        for trace in tracemap.contents:
            fibid = trace.fibid
            idx = trace.fibid - 1

            if trace.valid:
                self.logger.info('-' * 52)
                self.logger.info('Starting row %d, fibid %d', idx, fibid)

                row = rss[idx]

                fxpeaks, sxpeaks = find_fxpeaks(
                    sp=row,
                    times_sigma_threshold=0.0,
                    minimum_threshold=0,
                    nwinwidth_initial=7,
                    nwinwidth_refined=5,
                    npix_avoid_border=6,
                    nbrightlines=nlines,
                    sigma_gaussian_filtering=0,
                    minimum_gaussian_filtering=0
                )
                self.logger.info('number of peaks (expected): %s',
                                 str(nlines))
                self.logger.info('number of peaks (found)...: %d',
                                 len(fxpeaks))

                try:

                    # use channels (pixels from 1 to naxis1)
                    xchannel = fxpeaks + 1.0

                    list_of_wvfeatures = arccalibration_direct(
                        wv_master=wv_master,
                        ntriplets_master=ntriplets_master,
                        ratios_master_sorted=ratios_master_sorted,
                        triplets_master_sorted_list=triplets_master_sorted_list,
                        xpos_arc=xchannel,
                        naxis1_arc=naxis1,
                        crpix1=crpix1,
                        wv_ini_search=wv_ini_search,
                        wv_end_search=wv_end_search,
                        error_xpos_arc=3.0, # initially: 2.0
                        times_sigma_r=3.0,
                        frac_triplets_for_sum=0.50,
                        times_sigma_theil_sen=10.0,
                        poly_degree_wfit=poldeg_initial,
                        times_sigma_polfilt=10.0,
                        times_sigma_cook=10.0,
                        times_sigma_inclusion=10.0,
                        debugplot=debugplot
                    )

                    self.logger.info('Solution for row %d completed', idx)
                    self.logger.info('Fitting solution for row %d', idx)
                    solution_wv = fit_list_of_wvfeatures(
                            list_of_wvfeatures,
                            naxis1_arc=naxis1,
                            crpix1=crpix1,
                            poly_degree_wfit=poldeg_initial,
                            weighted=False,
                            debugplot=0,
                            plot_title=None
                        )

                    self.logger.info('linear crval1, cdelt1: %f %f',
                                     solution_wv.cr_linear.crval,
                                     solution_wv.cr_linear.cdelt)

                    self.logger.info('fitted coefficients %s',
                                     solution_wv.coeff)

                    # store results for plotting
                    plot_tracenumber.append(fibid)
                    plot_npeaksfound.append(len(fxpeaks))
                    plot_crval1.append(solution_wv.cr_linear.crval)
                    plot_cdelt1.append(solution_wv.cr_linear.cdelt)
                    plot_coeff.append(solution_wv.coeff)

                    trace_pol = trace.polynomial
                    # Update feature with measurements of Y coord in original
                    # image
                    # Peak and FWHM in RSS
                    for feature in solution_wv.features:
                        # Compute Y
                        feature.ypos = trace_pol(feature.xpos)
                        # FIXME: check here FITS vs PYTHON coordinates, etc
                        peak_int = int(feature.xpos)
                        try:
                            peak, fwhm = self.calc_fwhm_of_line(row, peak_int,
                                                                lwidth=20)
                        except Exception as error:
                            self.logger.warning("%s", error)
                            self.logger.warning('error in feature %s', feature)
                            # workaround
                            peak = row[peak_int]
                            fwhm = 0.0
                        # I would call this peak instead...
                        feature.peak = peak
                        feature.fwhm = fwhm

                    new = FiberSolutionArcCalibration(fibid, solution_wv)
                    initial_data_wlcalib.contents.append(new)

                except (ValueError, TypeError, IndexError) as error:
                    self.logger.warning("%s", error)
                    self.logger.warning('problem in row %d, fibid %d', idx, fibid)
                    initial_data_wlcalib.error_fitting.append(fibid)
                    error_contador += 1

            else:
                self.logger.info('skipping row %d, fibid %d, not extracted',
                                 idx, fibid)
                missing_fib += 1
                initial_data_wlcalib.missing_fibers.append(fibid)

        self.logger.info('Errors in fitting: %s', error_contador)
        self.logger.info('Missing fibers: %s', missing_fib)

        # save PDF file with plots in working directory
        if self.intermediate_results:
            from numina.array.display.matplotlib_qt import plt
            from matplotlib.backends.backend_pdf import PdfPages
            pdf = PdfPages('wavecal_iter1.pdf')
            for dumplot in zip([plot_npeaksfound, plot_crval1, plot_cdelt1],
                               ['number of peaks found',
                                'linear CRVAL1 ' + r'($\AA$)',
                                'linear_CDELT1 ' + r'($\AA$/pixel)']):
                ax = ximplotxy(plot_tracenumber, dumplot[0],
                               xlabel='fiber number', ylabel=dumplot[1],
                               linestyle='', marker='.', color='C0',
                               show=False)
                pdf.savefig()
                plt.close()
            for ideg in range(poldeg_initial + 1):
                dumplot = [coef[ideg] for coef in plot_coeff]
                ax = ximplotxy(plot_tracenumber, dumplot,
                               xlabel='fiber number',
                               ylabel='coef[' + str(ideg) + ']',
                               linestyle='', marker='.', color='C0',
                               show=False)
                pdf.savefig()
                plt.close()
            pdf.close()

        self.logger.info('Generating fwhm_image...')
        image = self.generate_fwhm_image(initial_data_wlcalib.contents)
        fwhm_image = fits.PrimaryHDU(image)
        fwhm_hdulist = fits.HDUList([fwhm_image])

        if refine_wv_calibration:
            self.logger.info('Improving wavelength calibration...')
            # model polynomial coefficients vs. fiber number using
            # previous results stored in data_wlcalib
            list_poly_vs_fiber = self.model_coeff_vs_fiber(
                initial_data_wlcalib, poldeg_initial,
                times_sigma_reject=5)
            # recompute data_wlcalib from scratch
            missing_fib = 0
            error_contador = 0
            data_wlcalib = WavelengthCalibration(instrument='MEGARA')
            data_wlcalib.total_fibers = tracemap.total_fibers
            plot_tracenumber = []
            plot_npointseff = []
            plot_residualstd = []
            plot_crval1 = []
            plot_cdelt1 = []
            plot_coeff = []

            # output PDF with refined fits
            if store_pdf_with_refined_fits == 1:
                if self.intermediate_results:
                    if not os.path.exists('refined_wavecal'):
                        try:
                            os.makedirs('refined_wavecal')
                        except OSError as exc:  # Guard against race condition
                            if exc.errno != errno.EEXIST:
                                raise
            # refine wavelength calibration polynomial for each valid fiber
            for trace in tracemap.contents:
                fibid = trace.fibid
                idx = trace.fibid - 1
                if trace.valid:
                    self.logger.info('-' * 52)
                    self.logger.info('Starting row %d, fibid %d', idx, fibid)
                    # select spectrum for current fiber
                    row = rss[idx]
                    naxis1 = row.shape[0]
                    # estimate polynomial coefficients for current fiber
                    coeff = numpy.zeros(poldeg_initial + 1)
                    for k in range(poldeg_initial + 1):
                        dumpol = list_poly_vs_fiber[k]
                        coeff[k] = dumpol(fibid)
                    wlpol = numpy.polynomial.Polynomial(coeff)
                    # refine polynomial fit using the full set of arc lines
                    # in master list (and save output PDF file when requested)
                    if store_pdf_with_refined_fits == 1:
                        if self.intermediate_results:
                            from matplotlib.backends.backend_pdf import PdfPages
                            plottitle = 'fiber #{0:03d}'.format(fibid)
                            pdf = PdfPages(
                                'refined_wavecal/{0:03d}.pdf'.format(fibid)
                            )
                        else:
                            plottitle = None
                            pdf = None
                    else:
                        plottitle = None
                        pdf = None
                    poly_refined, yres_summary  = \
                        refine_arccalibration(sp=row,
                                              poly_initial=wlpol,
                                              wv_master=wv_master_all,
                                              poldeg=poldeg_refined,
                                              plottitle=plottitle,
                                              ylogscale=True,
                                              pdf=pdf)
                    if pdf is not None:
                        from numina.array.display.matplotlib_qt import plt
                        plt.close()
                        pdf.close()
                    if poly_refined != numpy.polynomial.Polynomial([0.0]):
                        npoints_eff = yres_summary['npoints']
                        residual_std = yres_summary['robust_std']
                        # compute approximate linear values
                        crmin1_linear = poly_refined(1)
                        crmax1_linear = poly_refined(naxis1)
                        cdelt1_linear = (crmax1_linear - crmin1_linear) / \
                                        (naxis1 - 1)
                        self.logger.info('linear crval1, cdelt1: %f %f',
                                         crmin1_linear, cdelt1_linear)
                        self.logger.info('fitted coefficients %s',
                                         poly_refined.coef)
                        self.logger.info('npoints_eff, residual_std: %d %f',
                                         npoints_eff, residual_std)
                        cr_linear = CrLinear(
                            crpix=1.0,
                            crval=crmin1_linear,
                            crmin=crmin1_linear,
                            crmax=crmax1_linear,
                            cdelt=cdelt1_linear
                        )
                        solution_wv = SolutionArcCalibration(
                            features=[],  # empty list!
                            coeff=poly_refined.coef,
                            residual_std=residual_std,
                            cr_linear=cr_linear
                        )
                        solution_wv.npoints_eff = npoints_eff  # add also this
                        new = FiberSolutionArcCalibration(fibid, solution_wv)
                        data_wlcalib.contents.append(new)
                        # store results for plotting
                        plot_tracenumber.append(fibid)
                        plot_npointseff.append(npoints_eff)
                        plot_residualstd.append(residual_std)
                        plot_crval1.append(crmin1_linear)
                        plot_cdelt1.append(cdelt1_linear)
                        plot_coeff.append(poly_refined.coef)
                    else:
                        self.logger.error('error in row %d, fibid %d',
                                          idx, fibid)
                        data_wlcalib.error_fitting.append(fibid)
                else:
                    self.logger.info('skipping row %d, fibid %d, not extracted',
                                     idx, fibid)
                    missing_fib += 1
                    data_wlcalib.missing_fibers.append(fibid)

            self.logger.info('Errors in fitting: %s', error_contador)
            self.logger.info('Missing fibers: %s', missing_fib)

            # save PDF file with plots in working directory
            if self.intermediate_results:
                from numina.array.display.matplotlib_qt import plt
                from matplotlib.backends.backend_pdf import PdfPages
                pdf = PdfPages('wavecal_iter2.pdf')
                for dumplot in zip(
                        [plot_npointseff, plot_residualstd, plot_crval1,
                         plot_cdelt1],
                        ['effective number of lines found',
                         'residual std ' + r'($\AA$)',
                         'linear CRVAL1 ' + r'($\AA$)',
                         'linear_CDELT1 ' + r'($\AA$/pixel)']):
                    ax = ximplotxy(plot_tracenumber, dumplot[0],
                                   xlabel='fiber number', ylabel=dumplot[1],
                                   linestyle='', marker='.', color='C0',
                                   show=False)
                    pdf.savefig()
                    plt.close()
                for ideg in range(poldeg_refined + 1):
                    dumplot = [coef[ideg] for coef in plot_coeff]
                    ax = ximplotxy(plot_tracenumber, dumplot,
                                   xlabel='fiber number',
                                   ylabel='coef[' + str(ideg) + ']',
                                   linestyle='', marker='.', color='C0',
                                   show=False)
                    pdf.savefig()
                    plt.close()
                pdf.close()
        else:
            data_wlcalib = None

        self.logger.info('End arc calibration')

        return initial_data_wlcalib, data_wlcalib, fwhm_hdulist

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

    def model_coeff_vs_fiber(self, data_wlcalib, poldeg,
                             times_sigma_reject=5):
        """Model polynomial coefficients vs. fiber number.

        For each polynomial coefficient, a smooth polynomial dependence
        with fiber number is also computed, rejecting information from
        fibers which coefficients depart from that smooth variation.
        """

        if self.intermediate_results:
            from numina.array.display.matplotlib_qt import plt
            from matplotlib.backends.backend_pdf import PdfPages
            pdf = PdfPages('wavecal_refine_iter1.pdf')
            local_debugplot = 11
        else:
            pdf = None
            local_debugplot = 0

        list_fibid = []
        list_coeffs = []
        for item in (data_wlcalib.contents):
            list_fibid.append(item.fibid)
            if len(item.solution.coeff) != poldeg + 1:
                raise ValueError('Unexpected number of polynomial '
                                 'coefficients')
            list_coeffs.append(item.solution.coeff)

        # determine bad fits from each independent polynomial coefficient
        # (bad fits correspond to unexpected coefficient values for any of
        # the coefficients; i.e., the number of bad fits increases as we
        # examine different coefficients)
        poldeg_coeff_vs_fiber = 5
        reject_all = None  # avoid PyCharm warning
        fibid = numpy.array(list_fibid)
        for i in range(poldeg + 1):
            coeff = numpy.array([coeff[i] for coeff in list_coeffs])
            poly, yres, reject = polfit_residuals_with_sigma_rejection(
                x=fibid,
                y=coeff,
                deg=poldeg_coeff_vs_fiber,
                times_sigma_reject=times_sigma_reject,
            )
            if pdf is not None:
                polfit_residuals(
                    x=fibid,
                    y=coeff,
                    deg=poldeg_coeff_vs_fiber,
                    reject=reject,
                    xlabel='fibid',
                    ylabel='coeff a_' + str(i),
                    title='Identifying bad fits',
                    show=False,
                    debugplot=local_debugplot
                )
                pdf.savefig()
                plt.close()

            if i == 0:
                # initialize bad fits
                reject_all = numpy.copy(reject)
            else:
                # add new bad fits
                reject_all = numpy.logical_or(reject_all, reject)
            dumlabel = 'coeff a_' + str(i) + ': nreject=' + \
                       str(sum(reject_all))
            self.logger.info(dumlabel)
            self.logger.info(fibid[reject_all])

        # determine new fits excluding all fibers with bad fits
        list_poly_vs_fiber = []
        for i in range(poldeg + 1):
            coeff = numpy.array([coeff[i] for coeff in list_coeffs])
            poly, yres = polfit_residuals(
                x=fibid,
                y=coeff,
                deg=poldeg_coeff_vs_fiber,
                reject=reject_all,
                xlabel='fibid',
                ylabel='coeff a_' + str(i),
                title='Computing filtered fits',
                show=False,
                debugplot=local_debugplot
            )
            if pdf is not None:
                pdf.savefig()
                plt.close()
            list_poly_vs_fiber.append(poly)

        self.logger.info("list_poly_vs_fiber:")
        for i in range(poldeg + 1):
            self.logger.info(list_poly_vs_fiber[i])

        if pdf is not None:
            pdf.close()

        return list_poly_vs_fiber
