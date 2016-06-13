#
# Copyright 2011-2015 Universidad Complutense de Madrid
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
import matplotlib.pyplot as plt

from numina.array.trace.traces import trace
from numina.core import Product
from numina.core.requirements import ObservationResultRequirement

from megaradrp.products import MasterFiberFlatFrame, TraceMap
from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.requirements import MasterBiasRequirement, MasterBPMRequirement
from megaradrp.requirements import MasterDarkRequirement
from megaradrp.trace.traces import init_traces_ex

from skimage.filters import threshold_otsu

_logger = logging.getLogger('numina.recipes.megara')


class TraceMapRecipe(MegaraBaseRecipe):

    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_bpm = MasterBPMRequirement()

    fiberflat_frame = Product(MasterFiberFlatFrame)
    master_traces = Product(TraceMap)

    def __init__(self):
        super(TraceMapRecipe, self).__init__(version="0.1.0")

    def run(self, rinput):
        parameters = self.get_parameters(rinput)
        reduced = self.bias_process_common(rinput.obresult, parameters)

        data = reduced[0].data
        # For a given VPH, the position of the borders of the boxes
        # depend on position
        # For our current VPH
        current_vph = rinput.obresult.tags['vph']
        cstart = rinput.obresult.configuration.values['box']['boxcol']
        box_borders = rinput.obresult.configuration.values['box'][current_vph]

        hs = 3
        step1 = 2
        poldeg = 5
        maxdis1 = 2.0

        _logger.info('estimate background in column %i', cstart)
        background = estimate_background(data, center=cstart, hs=hs, boxref=box_borders)
        _logger.info('background level is %f', background)

        _logger.info('find peaks in column %i', cstart)

        central_peaks = init_traces_ex(data, center=cstart, hs=hs, box_borders=box_borders, tol=1.63)

        # The byteswapping is required by the cython module
        if data.dtype.byteorder != '=':
            _logger.debug('byteswapping image')
            image2 = data.byteswap().newbyteorder()
        else:
            image2 = data

        tracelist = []
        _logger.info('trace peaks')
        for dtrace in central_peaks.values():
            # FIXME, for traces, the background must be local, the background
            # in the center is not always good
            local_trace_background = 300 # background
            if dtrace.start:
                mm = trace(image2, x=cstart, y=dtrace.start[1], step=step1,
                         hs=hs, background=local_trace_background, maxdis=maxdis1)
                if False:
                    plt.plot(mm[:,0], mm[:,1])
                    plt.savefig('trace-xy-%d.png' % dtrace.fibid)
                    plt.close()
                    plt.plot(mm[:,0], mm[:,2])
                    plt.savefig('trace-xz-%d.png' % dtrace.fibid)
                    plt.close()
                if len(mm) < poldeg + 1:
                    _logger.warning('in fibid %d, only %d points to fit pol of degree %d',
                                    dtrace.fibid, len(mm), poldeg)
                    pfit = numpy.array([])
                else:
                    pfit = numpy.polyfit(mm[:,0], mm[:,1], deg=poldeg)

                start = mm[0, 0]
                stop = mm[-1,0]
            else:
                pfit = numpy.array([])
                start = cstart
                stop = cstart

            tracelist.append({'fibid': dtrace.fibid, 'boxid': dtrace.boxid,
                              'start':int(start), 'stop':int(stop),
                              'fitparms': pfit.tolist()})

        return self.create_result(fiberflat_frame=reduced,
                                  master_traces=tracelist)


def estimate_background(image, center, hs, boxref):
    """Estimate background from values in boxes between fibers"""

    cut_region = slice(center-hs, center+hs)
    cut = image[boxref, cut_region]

    colcut = cut.mean(axis=1)

    return threshold_otsu(colcut)
