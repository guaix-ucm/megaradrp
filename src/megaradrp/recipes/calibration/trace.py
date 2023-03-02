#
# Copyright 2011-2021 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Fiber tracing Recipe."""

from __future__ import division, print_function

import logging
import warnings

import numpy
import numpy.polynomial.polynomial as nppol
from numina.array.peaks.peakdet import refine_peaks
from numina.array.trace.traces import trace, tracing_limits
from numina.core import Result, Parameter
import matplotlib.pyplot as plt

import numina.types.qc as qc
from numina.array import combine
from numina.array.wavecalib.crosscorrelation import cosinebell
from numina.array.wavecalib.crosscorrelation import convolve_comb_lines
from skimage.filters import threshold_otsu
from skimage.feature import peak_local_max
from scipy.ndimage.filters import minimum_filter
from numina.frame.utils import copy_img

from megaradrp.processing.aperture import ApertureExtractor
from megaradrp.processing.combine import basic_processing_with_combination
from megaradrp.products import TraceMap
from megaradrp.products.tracemap import GeometricTrace
from megaradrp.ntypes import ProcessedImage, ProcessedRSS
from megaradrp.core.recipe import MegaraBaseRecipe
import megaradrp.requirements as reqs
import megaradrp.products
import megaradrp.processing.fibermatch as fibermatch
from megaradrp.instrument import vph_thr
from megaradrp.instrument.focalplane import FocalPlaneConf


class TraceMapRecipe(MegaraBaseRecipe):
    """Provides tracing information from continuum flat images.

    This recipe process a set of continuum flat images obtained in
    *Trace Map* mode and returns the tracing information required
    to perform fiber extraction in other recipes. The recipe also
    returns the result of processing the input images upto dark correction.

    See Also
    --------
    megaradrp.products.tracemap.TraceMap: description of TraceMap product
    megaradrp.recipes.calibration.modelmap.ModelMapRecipe: description of ModelMap recipe
    numina.array.trace.traces: tracing algorithm
    megaradrp.instrument.configs: instrument configuration

    Notes
    -----
    Images provided in `obresult` are trimmed and corrected from overscan,
    bad pixel mask (if `master_bpm` is not None), bias and dark current
    (if `master_dark` is not None).
    Images thus corrected are the stacked using the median.

    The result of the combination is saved as an intermediate result, named
    'reduced_image.fits'. This combined image is also returned in the field
    `reduced_image` of the recipe result and will be used for
    tracing the position of the fibers.

    The fibers are grouped in packs of different numbers of fibers. To match
    the traces in the image with the corresponding fibers is neccessary
    to know how fibers are packed and where the different groups of fibers
    appear in the detector. This information is provided by the fields
    'pseudoslit.boxes' and 'pseudoslit.boxes_positions' of the instrument
    configuration.

    Using the column reference provided by 'pseudoslit.boxes_positions', peaks
    are detected (using an average of 7 columns) and matched to the layout
    of fibers provided by 'pseudoslit.boxes_positions'. Fibers without a matching peak
    are counted and their ids stored in the final `master_traces` object.

    Once the peaks in the reference column are found, each one is traced
    until the border of the image is reached. The trace may be lost before
    reaching the border. In all cases, the beginning and the end of the trace
    are stored.

    The Y position of the trace is fitted to a polynomial
    of degree `polynomial_degree`. The coefficients of the polynomial are
    stored in the final `master_traces` object.

    """
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_bpm = reqs.MasterBPMRequirement()
    polynomial_degree = Parameter(5, 'Polynomial degree of trace fitting')
    relative_threshold = Parameter(0.3, 'Threshold for peak detection')
    debug_plot = Parameter(0, 'Save intermediate tracing plots')

    reduced_image = Result(ProcessedImage)
    reduced_rss = Result(ProcessedRSS)
    master_traces = Result(TraceMap)

    def run_qc(self, recipe_input, recipe_result):
        """Run quality control checks"""
        self.logger.info('start trace recipe QC')

        self.check_qc_tracemap(recipe_result.master_traces)

        recipe_result.qc = recipe_result.master_traces.quality_control
        self.logger.info('Result QC is %s', recipe_result.qc)
        self.logger.info('end trace recipe QC')
        return recipe_result

    def check_qc_tracemap(self, tracemap):
        # check full range in all fibers
        if tracemap.tags['insmode'] == 'LCB':
            nfibers = 623
            ignored_ids = [623]
        else:
            nfibers = 644
            ignored_ids = [635]

        cost_ids = nfibers * [0]

        start_min, end_max = tracemap.expected_range

        for trace in tracemap.contents:
            idx = trace.fibid - 1
            if trace.fibid in ignored_ids:
                continue
            if trace.start > start_min:
                cost_ids[idx] += 1
                msg = f'In fiber {trace.fibid}, trace start > {start_min}'
                warnings.warn(msg)
            if trace.stop < end_max:
                cost_ids[idx] += 1
                msg = f'In fiber {trace.fibid}, trace end < {end_max}'
                warnings.warn(msg)

        total = sum(cost_ids)

        if total > 300:
            tracemap.quality_control = qc.QC.BAD
        elif total > 0:
            tracemap.quality_control = qc.QC.PARTIAL
        else:
            tracemap.quality_control = qc.QC.GOOD

    # @numina.core.validator.validate
    def run(self, rinput):
        """Execute the recipe.

        Parameters
        ----------
        rinput : TraceMapRecipe.RecipeInput

        Returns
        -------
        TraceMapRecipe.RecipeResult

        Raises
        ------
        ValidationError
              if the recipe input is invalid

        """
        self.logger.info('start trace spectra recipe')

        obresult = rinput.obresult

        obresult_meta = obresult.metadata_with(self.datamodel)

        debug_plot = rinput.debug_plot if self.intermediate_results else 0

        self.logger.info('start basic reduction')
        flow = self.init_filters(rinput, obresult.configuration)
        reduced = basic_processing_with_combination(rinput, flow, method=combine.median)
        self.logger.info('end basic reduction')

        self.save_intermediate_img(reduced, 'reduced_image.fits')

        # insconf = obresult.configuration
        insconf = obresult.profile

        boxes = insconf.get_property('pseudoslit.boxes')
        values = insconf.get_property('pseudoslit.boxes_positions')
        cstart0 = values['ref_column']
        box_borders0 = values['positions']

        box_borders, cstart = self.refine_boxes_from_image(reduced, box_borders0, cstart0)

        self.logger.debug("original boxes: %s", box_borders0)
        self.logger.debug("refined boxes: %s", box_borders)

        current_vph = reduced[0].header['vph']
        current_insmode = reduced[0].header['insmode']

        if current_insmode in vph_thr and current_vph in vph_thr[current_insmode]:
            threshold = vph_thr[current_insmode][current_vph]
            self.logger.info('rel threshold for %s is %4.2f', current_vph, threshold)
        else:
            threshold = rinput.relative_threshold
            self.logger.info('rel threshold not defined for %s, using %4.2f', current_vph, threshold)

        final = megaradrp.products.TraceMap(instrument=obresult.instrument)
        fp_conf = FocalPlaneConf.from_img(reduced)
        # As of 2019-10-15, headers do not contain information
        # about inactive fibers, i.e. all have active=True
        # inactive_fibers = fp_conf.inactive_fibers()
        # We do it manually
        if fp_conf.name == 'LCB':
            inactive_fibers = [623]
        else:
            inactive_fibers = [635]

        final.total_fibers = fp_conf.nfibers
        final.tags = self.extract_tags_from_ref(reduced, final.tag_names(), base=obresult.labels)
        final.boxes_positions = box_borders
        final.ref_column = cstart

        # Searching for peaks
        # step of
        step = 2
        # number of the columns to add
        hs = 3
        # Expected range of computed traces
        xx_start, xx_end = tracing_limits(reduced[0].shape[1], cstart, step, hs)
        final.expected_range = [xx_start, xx_end]

        final.update_metadata(self)
        final.update_metadata_origin(obresult_meta)
        # Temperature in Celsius with 2 decimals
        final.tags['temp'] = round(obresult_meta['info'][0]['temp'] - 273.15, 2)

        contents, error_fitting, missing_fibers = self.search_traces(
            reduced,
            boxes,
            box_borders,
            inactive_fibers=inactive_fibers,
            cstart=cstart,
            step=step,
            hs=hs,
            threshold=threshold,
            poldeg=rinput.polynomial_degree,
            debug_plot=debug_plot
        )

        final.contents = contents
        final.error_fitting = error_fitting
        final.missing_fibers = missing_fibers
        # Perform extraction with own traces
        calibrator_aper = ApertureExtractor(final, self.datamodel)
        reduced_copy = copy_img(reduced)
        reduced_rss = calibrator_aper(reduced_copy)

        if self.intermediate_results:
            with open('ds9.reg', 'w') as ds9reg:
                final.to_ds9_reg(ds9reg, rawimage=False,
                                 numpix=100, fibid_at=2048)

            with open('ds9_raw.reg', 'w') as ds9reg:
                final.to_ds9_reg(ds9reg, rawimage=True,
                                 numpix=100, fibid_at=2048)

        self.logger.info('end trace spectra recipe')
        return self.create_result(reduced_image=reduced,
                                  reduced_rss=reduced_rss,
                                  master_traces=final)

    def obtain_boxes_from_image(self, reduced, expected, npeaks, cstart=2000):
        from numina.array.peaks.peakdet import find_peaks_indexes
        col = cstart
        data = reduced[0].data
        rr = data[:, col-1:col+1].mean(axis=1)
        # standarize
        rr -= numpy.median(rr)
        rr /= rr.max()
        rr *= -1

        cb = cosinebell(len(rr), 0.10)
        cbr = cb * rr
        plt.plot(cbr)
        plt.show()
        xv = numpy.fft.fftfreq(len(cbr))
        yv = numpy.fft.fft(cbr)
        plt.xlim([0.0, 0.5])
        plt.semilogy(xv, numpy.abs(yv.real))
        plt.show()

        cut = abs(xv) > 0.1
        yv[cut] = 0
        res = numpy.fft.ifft(yv)
        final = res.real
        plt.plot(final)
        # trend = detrend(final)
        # plt.plot(final - trend)
        plt.show()

        idx = find_peaks_indexes(final, window_width=3, threshold=0.3, fpeak=1)

        # We expect  differentnumber of pekas in LCB/MOS
        # order by intensity
        peak_flux = final[idx]
        # Number of peaks must be >=18
        npeaks = npeaks + 1
        fidx = numpy.argsort(peak_flux)[:-(npeaks+1):-1]
        nidx = idx[fidx]
        nidxs = numpy.sort(nidx)

        plt.plot(final)
        # plt.scatter(idx, [0.9 for m in idx])
        plt.scatter(nidx, [0.95 for m in nidx], c='r')
        plt.scatter(expected, [1.0 for m in expected])
        plt.show()

        plt.scatter(expected, nidxs - expected)
        plt.show()
        return nidxs, col

    def refine_boxes_from_image(self, reduced, expected, cstart=2000, nsearch=20):
        """Refine boxes using a filtered Fourier image"""

        hs = 3
        # Cut freq in Fourier space
        cut_frec = 0.10
        # Cosine bell
        cos_cut = 0.10

        data = reduced[0].data
        rr = data[:, cstart-hs:cstart+hs].mean(axis=1)

        # standarize and Y flip
        rr -= numpy.median(rr)
        rr /= rr.max()
        rr *= -1

        cb = cosinebell(len(rr), cos_cut)
        cbr = cb * rr

        xv = numpy.fft.fftfreq(len(cbr))
        yv = numpy.fft.fft(cbr)

        # Filter freqs in Fourier space
        yv[abs(xv) > cut_frec] = 0

        res = numpy.fft.ifft(yv)
        final = res.real

        # initial determination of global offset (integer number) by
        # cross-correlating an artificial spectrum (with lines placed at
        # the expected locations of the frontiers)
        expected_arr = numpy.array(expected)
        comb_lines_flux = numpy.ones_like(expected_arr)
        naxis1 = len(final)
        ioffset_max = 100  # pixels
        xcorr = []
        ycorr = []
        xwave = None   # avoid PyCharm warning
        sp_comb_lines0 = None   # avoid PyCharm warning
        for ioffset in range(-ioffset_max, ioffset_max + 1):
            xwave, sp_comb_lines = convolve_comb_lines(
                lines_wave=expected_arr + ioffset,
                lines_flux=comb_lines_flux,
                sigma=3,
                crpix1=0, crval1=0, cdelt1=1, naxis1=naxis1
            )
            factor = max(final) / max(sp_comb_lines)
            sp_comb_lines = sp_comb_lines * factor
            fpeak = numpy.sum(sp_comb_lines * final)
            xcorr.append(ioffset)
            ycorr.append(fpeak)
            if ioffset == 0:
                sp_comb_lines0 = sp_comb_lines.copy()  # store for plot

        # initial offset
        ioffset = xcorr[ycorr.index(max(ycorr))]

        # auxiliary plot showing the cross-correlation work
        if self.intermediate_results:
            fig, ax = plt.subplots(ncols=1, nrows=1)
            ax.plot(xcorr, ycorr, 'o-')
            ax.set_xlabel('ioffset')
            ax.set_ylabel('peak of corrrelation function')
            ax.axvline(ioffset, linestyle=':', color='C1')
            ax.set_title(f'optimal initial offset: {ioffset}')
            plt.savefig('offset_borders_cross.png')
            plt.close()

        # using the initial offset, refine the peak search around the new
        # expected location, looking for a maximum in +/- nsearch pixels
        refined = expected[:]    # initialize list with the same length
        for ibox, box in enumerate(expected):
            box_ini = box - nsearch + ioffset
            box_end = box + nsearch + 1 + ioffset
            iargmax = final[box_ini:box_end].argmax()
            refined[ibox] = iargmax + box - nsearch + ioffset

        # auxiliary plot showing the initial and final frontier locations
        if self.intermediate_results:
            fig, ax = plt.subplots(ncols=1, nrows=1)
            ax.plot(final, label=f'cross section at x={cstart}')
            ax.plot(xwave, sp_comb_lines0, label='expected location of frontiers')
            for idum, item in enumerate(expected):
                if idum == 0:
                    label = 'refined frontier location'
                else:
                    label = None
                ax.plot(refined, final[refined], 'ro', label=label)
            ax.legend()
            ax.set_title(f'optimal initial offset: {ioffset}')
            ax.set_xlabel('pixel number - 1')
            ax.set_ylabel('inverted normalized signal')
            plt.savefig('frontiers_between_pseudoslits.png')
            plt.close()

        return refined, cstart

    def search_traces(self, reduced, boxes, box_borders, inactive_fibers=None, cstart=2000,
                      threshold=0.3, poldeg=5, step=2, hs=3, debug_plot=0):

        data = reduced[0].data
        if inactive_fibers is None:
            inactive_fibers = []
        tol = 1.63

        self.logger.info('search for traces')

        self.logger.info('estimate background in ref column %i', cstart)
        background = estimate_background(data, center=cstart, hs=hs, boxref=box_borders)
        self.logger.info('background level is %f', background)

        self.logger.info('find peaks in reference column %i', cstart)

        central_peaks = init_traces(
            data,
            center=cstart,
            hs=hs,
            boxes=boxes,
            box_borders=box_borders,
            tol=tol,
            threshold=threshold,
            debug_plot=debug_plot
        )

        # The byteswapping is required by the cython module
        if data.dtype.byteorder != '=':
            self.logger.debug('byteswapping image')
            image2 = data.byteswap().newbyteorder()
        else:
            image2 = data

        maxdis = 2.0

        contents = []
        error_fitting = []
        missing_fibers = []
        self.logger.info('trace peaks from references')
        for dtrace in central_peaks:
            # FIXME, for traces, the background must be local
            # the background in the center is not always good
            local_trace_background = 300  # background

            self.logger.debug('trace fiber %d', dtrace.fibid)
            conf_ok = dtrace.fibid not in inactive_fibers
            peak_ok = dtrace.start is not None

            pfit = numpy.array([])
            start = cstart
            stop = cstart

            if peak_ok:
                if not conf_ok:
                    error_fitting.append(dtrace.fibid)
                    self.logger.warning('found fibid %d, expected to be missing', dtrace.fibid)
                else:

                    mm = trace(image2, x=cstart, y=dtrace.start[1], step=step,
                               hs=hs, background=local_trace_background, maxdis=maxdis)

                    if debug_plot:
                        plt.plot(mm[:, 0], mm[:, 1], '.')
                        plt.savefig(f'trace-xy-{dtrace.fibid:03d}.png')
                        plt.close()
                        plt.plot(mm[:, 0], mm[:, 2], '.')
                        plt.savefig(f'trace-xz-{dtrace.fibid:03d}.png')
                        plt.close()
                    if len(mm) < poldeg + 1:
                        self.logger.warning('in fibid %d, only %d points to fit pol of degree %d',
                                            dtrace.fibid, len(mm), poldeg)
                        pfit = numpy.array([])
                    else:
                        pfit = nppol.polyfit(mm[:, 0], mm[:, 1], deg=poldeg)

                    start = mm[0, 0]
                    stop = mm[-1, 0]
                    self.logger.debug('trace start %d  stop %d', int(start), int(stop))
            else:
                if conf_ok:
                    self.logger.warning('error tracing fibid %d', dtrace.fibid)
                    error_fitting.append(dtrace.fibid)
                else:
                    self.logger.debug('expected missing fibid %d', dtrace.fibid)
                    missing_fibers.append(dtrace.fibid)

            this_trace = GeometricTrace(
                fibid=dtrace.fibid,
                boxid=dtrace.boxid,
                start=int(start),
                stop=int(stop),
                fitparms=pfit.tolist()
            )
            contents.append(this_trace)

        return contents, error_fitting, missing_fibers


def estimate_background(image, center, hs, boxref):
    """Estimate background from values in boxes between fibers"""

    cut_region = slice(center-hs, center+hs)
    cut = image[boxref, cut_region]

    colcut = cut.mean(axis=1)

    return threshold_otsu(colcut)


class FiberTraceInfo(object):
    def __init__(self, fibid, boxid):
        self.boxid = boxid
        self.fibid = fibid
        self.start = None


def init_traces(image, center, hs, boxes, box_borders, tol=1.5, threshold=0.37, debug_plot=0):

    _logger = logging.getLogger(__name__)

    cut_region = slice(center-hs, center+hs)
    cut = image[:, cut_region]
    colcut = cut.mean(axis=1)
    # trace local background with a min filter
    mincol = minimum_filter(colcut, size=7)

    _logger.debug('initial pairing fibers in column %d', center)

    fiber_traces = []
    total_peaks = 0
    lastid = 0
    counted_fibers = 0
    boxes_with_missing_fibers = []

    ioffset = 0  # correction to be applied to each successive fiber block

    for boxid, box in enumerate(boxes):
        nfibers = box['nfibers']
        mfibers = box.get('missing', [])
        nfibers_max = nfibers - len(mfibers)
        sfibers = box.get('skipcount', [])
        _logger.debug('pseudoslit box: %s, id: %d', box['name'], boxid)
        _logger.debug('nfibers: %d, missing: %s', nfibers, mfibers)

        counted_fibers += nfibers
        b1 = int(box_borders[boxid])
        b2 = int(box_borders[boxid + 1])
        _logger.debug('initial box borders: %s %s', b1, b2)
        _logger.debug('offset for borders: %s', ioffset)
        b1 += ioffset
        b2 += ioffset
        _logger.debug('updated box borders: %s %s', b1, b2)
        borders = [b1, b2]

        region = colcut[borders[0]:borders[1]+1]
        # Region for background computation
        # Remove the space of 1.5 fibers
        bb1 = b1 + 9
        bb2 = b2 - 9
        # a conservative background
        lowt = numpy.min(mincol[bb1:bb2+1])
        _logger.debug('conservative background in box is %f', lowt)

        # Find exactly the number of peaks expected
        _logger.debug('find %d peaks (max)', nfibers_max)
        ipeaks_int = peak_local_max(region, min_distance=3,
                                    threshold_abs=lowt,
                                    threshold_rel=threshold,
                                    num_peaks=nfibers_max,
                                    )[:, 0]

        npeaks = len(ipeaks_int)
        total_peaks += npeaks

        if npeaks == 0:
            # skip everything, go to next box
            boxes_with_missing_fibers.append((box['name'], nfibers))
            _logger.debug('no peaks detected, go to next box')
            continue

        # We always want the result sorted. The order changes in different versions
        # of scikit-image
        ipeaks_int.sort()
        ipeaks_float = refine_peaks(region, ipeaks_int, 3)[0]
        peaks_dist = numpy.diff(ipeaks_float)
        measured_scale = numpy.median(peaks_dist)
        _logger.debug('median distance between peaks is %s', measured_scale)

        peaks_y = numpy.ones((ipeaks_int.shape[0], 3))
        peaks_y[:, 0] = ipeaks_int + b1
        peaks_y[:, 1] = ipeaks_float + b1
        peaks_y[:, 2] = region[ipeaks_int]
        peak1 = peaks_y[0, 0]
        peak2 = peaks_y[-1, 0]
        border1 = peak1 - b1  # distance from left border to first peak
        border2 = b2 - peak2  # distance from last peak to right border
        delta_ioffset = int((border1 - border2) / 2)
        ioffset += delta_ioffset  # new offset to recenter next fiber block

        if debug_plot:
            plt.plot(numpy.arange(len(region))+b1, region)
            plt.plot(ipeaks_int+b1, region[ipeaks_int], 'r*')
            plt.xlabel('pixel number - 1')
            plt.ylabel('number of counts')
            plt.title(f'pseudoslit box: {box["name"]}, id: {boxid}')
            plt.savefig(f'central_cut_{boxid:02d}.png')
            plt.close()

        startid = lastid + 1
        fiber_model = fibermatch.generate_box_model(nfibers,
                                                    start=startid,
                                                    missing_relids=mfibers,
                                                    skip_fibids=sfibers
                                                    )
        lastid = fiber_model[-1].fibid
        _logger.debug('fibids %s - %s', startid, lastid)

        matched_peaks = fibermatch.count_peaks(peaks_y[:, 1])
        nmatched = len(matched_peaks)
        missing = nfibers - nmatched

        if missing != 0:
            boxes_with_missing_fibers.append((box['name'], missing))

        _logger.debug('matched %s, missing: %s', nmatched, missing)
        pos_solutions = fibermatch.complete_solutions(fiber_model, matched_peaks,
                                                      borders, scale=measured_scale)

        for fibid, match in fibermatch.iter_best_solution(fiber_model,
                                                          matched_peaks,
                                                          pos_solutions):
            fti = FiberTraceInfo(fibid, boxid)
            if match is not None:
                fti.start = (center, peaks_y[match, 1], peaks_y[match, 2])
            fiber_traces.append(fti)

    _logger.debug('total found peaks: %d', total_peaks)
    _logger.debug('total found + recovered peaks: %d', counted_fibers)
    if boxes_with_missing_fibers:
        for m1, n2 in boxes_with_missing_fibers:
            _logger.debug('missing %d fibers in box %s', n2, m1)
    return fiber_traces
