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

"""Fiber tracing Recipe."""

from __future__ import division, print_function

import logging

import numpy
import numpy.polynomial.polynomial as nppol
from numina.array.peaks.peakdet import refine_peaks
from numina.array.trace.traces import trace
from numina.core import Product, Parameter
from numina.core.requirements import ObservationResultRequirement
from numina.array import combine
import numina.core.validator
from skimage.filters import threshold_otsu
from skimage.feature import peak_local_max

from megaradrp.instrument.loader import build_instrument_config
from megaradrp.processing.combine import basic_processing_with_combination
from megaradrp.products import TraceMap
from megaradrp.products.tracemap import GeometricTrace
from megaradrp.types import ProcessedFrame
from megaradrp.core.recipe import MegaraBaseRecipe
import megaradrp.requirements as reqs
import megaradrp.products

from megaradrp.instrument import vph_thr


_logger = logging.getLogger(__name__)


class TraceMapRecipe(MegaraBaseRecipe):

    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_bpm = reqs.MasterBPMRequirement()
    polynomial_degree = Parameter(5, 'Polynomial degree of trace fitting')
    relative_threshold = Parameter(0.3, 'Threshold for peak detection')

    fiberflat_frame = Product(ProcessedFrame)
    master_traces = Product(TraceMap)

    def __init__(self):
        super(TraceMapRecipe, self).__init__(version="0.1.0")

    @numina.core.validator.validate
    def run(self, rinput):
        self.logger.info('start trace spectra recipe')

        obresult = rinput.obresult
        current_vph = obresult.tags['vph']

        self.logger.info('start basic reduction')
        flow = self.init_filters(rinput, obresult.configuration)
        reduced = basic_processing_with_combination(rinput, flow, method=combine.median)
        self.logger.info('end basic reduction')

        #insconf_key = "4fd05b24-2ed9-457b-b563-a3c618bb1d4c"
        #print(obresult.configuration)
        #insconf = build_instrument_config(insconf_key)
        #print(insconf)
        insconf = obresult.configuration
        values = insconf.get('pseudoslit.boxes_positions', **obresult.tags)
        cstart = values['ref_column']
        box_borders = values['positions']

        boxes = insconf.get('pseudoslit.boxes', **obresult.tags)

        if current_vph in vph_thr:
            threshold = vph_thr[current_vph]
            self.logger.info('rel threshold for %s is %4.2f', current_vph, threshold)
        else:
            threshold = rinput.relative_threshold
            self.logger.info('rel threshold not defined for %s, using %4.2f', current_vph, threshold)

        final = megaradrp.products.TraceMap(instrument=obresult.instrument)
        final.tags = obresult.tags

        contents, error_fitting = self.search_traces(
            reduced,
            boxes,
            box_borders,
            cstart=cstart,
            threshold=threshold,
            poldeg=rinput.polynomial_degree
        )

        final.contents = contents
        final.error_fitting = error_fitting
        self.logger.info('end trace spectra recipe')
        return self.create_result(fiberflat_frame=reduced,
                                  master_traces=final)

    def search_traces(self, reduced, boxes, box_borders, cstart=2000, threshold=0.3, poldeg=5, step=2):

        data = reduced[0].data

        hs = 3
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
            threshold=threshold
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
        self.logger.info('trace peaks from references')
        for dtrace in central_peaks:
            # FIXME, for traces, the background must be local
            # the background in the center is not always good
            local_trace_background = 300  # background

            self.logger.debug('trace fiber %d', dtrace.fibid)
            if dtrace.start:
                mm = trace(image2, x=cstart, y=dtrace.start[1], step=step,
                           hs=hs, background=local_trace_background, maxdis=maxdis)
                if False:
                    import matplotlib.pyplot as plt
                    plt.plot(mm[:, 0], mm[:, 1])
                    plt.savefig('trace-xy-%d.png' % dtrace.fibid)
                    plt.close()
                    plt.plot(mm[:, 0], mm[:, 2])
                    plt.savefig('trace-xz-%d.png' % dtrace.fibid)
                    plt.close()
                if len(mm) < poldeg + 1:
                    self.logger.warning('in fibid %d, only %d points to fit pol of degree %d',
                                        dtrace.fibid, len(mm), poldeg)
                    pfit = numpy.array([])
                else:
                    pfit = nppol.polyfit(mm[:, 0], mm[:, 1], deg=poldeg)

                start = mm[0, 0]
                stop = mm[-1, 0]
            else:
                pfit = numpy.array([])
                start = cstart
                stop = cstart
                error_fitting.append(dtrace.fibid)

            self.logger.debug('trace start %d  stop %d', int(start), int(stop))

            this_trace = GeometricTrace(
                fibid=dtrace.fibid,
                boxid=dtrace.boxid,
                start=int(start),
                stop=int(stop),
                fitparms=pfit.tolist()
            )
            contents.append(this_trace)

        return contents, error_fitting


def estimate_background(image, center, hs, boxref):
    """Estimate background from values in boxes between fibers"""

    cut_region = slice(center-hs, center+hs)
    cut = image[boxref, cut_region]

    colcut = cut.mean(axis=1)

    return threshold_otsu(colcut)


# FIXME: need a better place for this
# Moved from megaradrp.trace


class FiberTraceInfo(object):
    def __init__(self, fibid, boxid):
        self.boxid = boxid
        self.fibid = fibid
        self.start = None


def init_traces(image, center, hs, boxes, box_borders, tol=1.5, threshold=0.37):

    cut_region = slice(center-hs, center+hs)
    cut = image[:,cut_region]
    colcut = cut.mean(axis=1)

    counted_fibers = 0
    fiber_traces = []
    total_peaks = 0
    total_peaks_pos = []

    # ipeaks_int = peak_local_max(colcut, min_distance=2, threshold_rel=0.2)[:, 0]
    ipeaks_int = peak_local_max(colcut, min_distance=3, threshold_rel=threshold)[:, 0] # All VPH

    if False:
        import matplotlib.pyplot as plt
        plt.plot(colcut)
        plt.plot(ipeaks_int, colcut[ipeaks_int], 'r*')
        for border in box_borders:
            plt.axvline(border, color='k')
        plt.show()

    ipeaks_float = refine_peaks(colcut, ipeaks_int, 3)[0]
    peaks_y = numpy.ones((ipeaks_int.shape[0], 3))
    peaks_y[:, 0] = ipeaks_int
    peaks_y[:, 1] = ipeaks_float
    peaks_y[:, 2] = colcut[ipeaks_int]
    box_match = numpy.digitize(peaks_y[:, 0], box_borders)

    _logger.debug('initial pairing fibers in column %d', center)
    for boxid, box in enumerate(boxes):
        nfibers = box['nfibers']
        dist_b_fibs = (box_borders[boxid + 1] - box_borders[boxid]) / (nfibers + 2.0)
        mask_fibers = (box_match == (boxid + 1))
        # Peaks in this box
        thispeaks = peaks_y[mask_fibers]
        npeaks = len(thispeaks)
        total_peaks += npeaks
        for elem in thispeaks:
            total_peaks_pos.append(elem.tolist())

        _logger.debug('pseudoslit box: %s, id: %d', box['name'], boxid)
        # Start by matching the first peak
        # with the first fiber
        fid = 0
        current_peak = 0
        pairs_1 = [(fid, current_peak)]
        fid += 1

        scale = 1
        while (current_peak < npeaks - 1) and (fid < nfibers):
            # Expected distance to next fiber
            expected_distance = scale * dist_b_fibs
            # _logger.debug('expected %s current peak %s', expected_distance, current_peak)
            for idx in range(current_peak + 1, npeaks):
                distance = abs(thispeaks[idx, 1] - thispeaks[current_peak, 1])
                if abs(distance - expected_distance) <= tol:
                    # We have a match
                    # We could update
                    # dist_b_fibs = distance / scale
                    # But is not clear this is better

                    # Store this match
                    pairs_1.append((fid, idx))
                    current_peak = idx
                    # Next
                    scale = 1
                    break
            else:
                # This fiber has no match
                pairs_1.append((fid, None))
                # Try a fiber further away
                scale += 1
            # Match next fiber
            fid += 1

        _logger.debug('matched %s \t missing: %s', len(pairs_1),nfibers-len(pairs_1))
        remainig = nfibers - len(pairs_1)
        if remainig > 0:
            _logger.debug('we have to pair %d miising fibers', remainig)
            # Position of first match fiber

            # Position of last match fiber
            for fid, peakid in reversed(pairs_1):
                if peakid is not None:
                    last_matched_peak = peakid
                    last_matched_fiber = fid
                    break
            else:
                raise ValueError('None matched')
            _logger.debug('peaks: %s \t %s', thispeaks[0, 1], thispeaks[last_matched_peak, 1])
            _logger.debug('borders: %s \t %s', box_borders[boxid], box_borders[boxid+1])
            ldist = thispeaks[0, 1] - box_borders[boxid]
            rdist = box_borders[boxid + 1] - thispeaks[last_matched_peak, 1]
            lcap = ldist / dist_b_fibs - 1
            rcap = rdist / dist_b_fibs - 1
            _logger.debug('L distance %s \t %s', ldist, lcap)
            _logger.debug('R distance %s \t %s', rdist, rcap)
            lcapi = int(lcap + 0.5)
            rcapi = int(rcap + 0.5)

            on_r = rcapi <= lcapi
            mincap = min(lcapi, rcapi)
            maxcap = max(lcapi, rcapi)

            cap1 = min(mincap, remainig)
            cap2 = min(maxcap, remainig - cap1)
            cap3 = remainig - cap1 - cap2

            if cap3 > 0:
                _logger.debug('we dont have space %s fibers no allocated', cap3)

            if on_r:
                # Fill rcap fibers, then lcap
                capr = cap1
                capl = cap2
            else:
                capr = cap2
                capl = cap1

            addl = [(x, None) for x in range(-capl, 0)]
            addr = [(x, None) for x in range(last_matched_fiber + 1, last_matched_fiber + 1 + capr)]
            _logger.debug('add %s fibers on the right', capr)
            _logger.debug('add %s fibers on the left', capl)
            pairs_1 = addl + pairs_1 + addr

        for fibid, (relfibid, match) in enumerate(pairs_1, counted_fibers):
            fti = FiberTraceInfo(fibid + 1, boxid)
            if match is not None:
                fti.start = (center, thispeaks[match, 1], thispeaks[match, 2])
            fiber_traces.append(fti)
            # else:
            #     fiber_traces[fibid].start = (center, 0, 0)
        counted_fibers += nfibers

        # import matplotlib.pyplot as plt
        # plt.xlim([box_borders[boxid], box_borders[boxid + 1]])
        # plt.plot(colcut, 'b-')
        # plt.plot(thispeaks[:, 1], thispeaks[:, 2], 'ro')
        # plt.plot(peaks_y[:,1], peaks_y[:, 2], 'ro')
        # plt.title('Box %s' %box['id'])
        # plt.show()

    # import matplotlib.pyplot as plt
    # total_peaks_pos = np.array(total_peaks_pos)
    # plt.plot(colcut, 'b-')
    # plt.plot(total_peaks_pos[:, 1], total_peaks_pos[:, 2], 'ro')
    # plt.show()

    _logger.debug('total found peaks: %d', total_peaks)
    _logger.debug('total found + recovered peaks: %d', counted_fibers)

    return fiber_traces
