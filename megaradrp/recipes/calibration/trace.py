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

from astropy.io import fits

from numina.array.trace.traces import trace
from numina.core import Product
from numina.core.requirements import ObservationResultRequirement

from megaradrp.products import MasterFiberFlat, TraceMap
from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.requirements import MasterBiasRequirement, MasterBPMRequirement
from megaradrp.requirements import MasterDarkRequirement
from megaradrp.trace.traces import init_traces

_logger = logging.getLogger('numina.recipes.megara')


class TraceMapRecipe(MegaraBaseRecipe):

    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_bpm = MasterBPMRequirement()
    fiberflat_frame = Product(MasterFiberFlat)
    master_traces = Product(TraceMap)

    def __init__(self):
        super(TraceMapRecipe, self).__init__(
            version="0.1.0"
        )

    def run(self, rinput):

        parameters = self.get_parameters(rinput)
        reduced = self.bias_process_common(rinput.obresult, parameters)

        data = reduced[0].data
        cstart = 2000
        hs = 3
        step1 = 2
        background1 = 150.0
        maxdis1 = 2.0

        _logger.info('find peaks in column %i', cstart)

        central_peaks = init_traces(data, center=cstart, hs=hs,
                                background=background1)

        _logger.info(' %i peaks found', len(central_peaks))


        # The byteswapping is required by the cython module
        if data.dtype.byteorder != '=':
            _logger.debug('byteswapping image')
            image2 = data.byteswap().newbyteorder()
        else:
            image2 = data

        tracelist = []
        _logger.info('trace peaks')
        for dtrace in central_peaks.values():

            mm = trace(image2, x=cstart, y=dtrace.start[1], step=step1,
                         hs=hs, background=background1, maxdis=maxdis1)

            pfit = numpy.polyfit(mm[:,0], mm[:,1], deg=5)

            tracelist.append({'fibid': dtrace.fibid, 'boxid': dtrace.boxid,
                              'start':0, 'stop':4095,
                              'fitparms': pfit.tolist()})

        return self.create_result(fiberflat_frame=reduced,
                                  master_traces=tracelist)
