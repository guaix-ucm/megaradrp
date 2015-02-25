#
# Copyright 2011-2014 Universidad Complutense de Madrid
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

'''Calibration Recipes for Megara'''

from __future__ import division, print_function

import logging

import numpy
from astropy.io import fits

from numina.core import Product
from numina.core.products import ArrayType
from numina.core.requirements import ObservationResultRequirement
from numina.array.combine import median as c_median
from numina.flow import SerialFlow
from numina.flow.processing import BiasCorrector

from megaradrp.core import MegaraBaseRecipe
from megaradrp.core import OverscanCorrector, TrimImage
from megaradrp.core import peakdet
# from numina.logger import log_to_history

from megaradrp.products import MasterFiberFlat
from megaradrp.products import TraceMapType, TraceMap
from megaradrp.requirements import MasterBiasRequirement

from megaradrp.trace.traces import init_traces
from megaradrp.trace._traces import tracing  # @UnresolvedImport





import numpy as np
import math



        
_logger = logging.getLogger('numina.recipes.megara')


def other(a, b):
    start = wcs_to_pix(a)
    end = wcs_to_pix(b)

    next = start + 1
    prev = end - 1
    print(start, start-0.5, start+0.5, min(start+0.5, b)-a)
    for i in range(next, prev+1):
        print(i,  1.0)
    if end > start:
        print(end, end-0.5, end+0.5, b - (end-0.5))

def fill_other(data, a, b):
    start = wcs_to_pix(a)
    end = wcs_to_pix(b)
    data[start] = min(start+0.5, b)-a
    data[start+1:end] = 1.0
    if end > start:
        data[end] = b - (end-0.5)
    return data

def wcs_to_pix(x):
    return int(math.floor(x + 0.5))


def process_common(recipe, obresult, master_bias):
    _logger.info('starting prereduction')

    o_c = OverscanCorrector()
    t_i = TrimImage()

    with master_bias.open() as hdul:
        mbias = hdul[0].data.copy()
        b_c = BiasCorrector(mbias)

    basicflow = SerialFlow([o_c, t_i, b_c])

    cdata = []

    try:
        for frame in obresult.frames:
            hdulist = frame.open()
            hdulist = basicflow(hdulist)
            cdata.append(hdulist)

        _logger.info('stacking %d images using median', len(cdata))

        data = c_median([d[0].data for d in cdata], dtype='float32')
        template_header = cdata[0][0].header
        hdu = fits.PrimaryHDU(data[0], header=template_header)
    finally:
        for hdulist in cdata:
            hdulist.close()

    hdr = hdu.header
    hdr['IMGTYP'] = ('FIBER_FLAT', 'Image type')
    hdr['NUMTYP'] = ('MASTER_FIBER_FLAT', 'Data product type')
    hdr = recipe.set_base_headers(hdr)
    hdr['CCDMEAN'] = data[0].mean()

    varhdu = fits.ImageHDU(data[1], name='VARIANCE')
    num = fits.ImageHDU(data[2], name='MAP')
    result = fits.HDUList([hdu, varhdu, num])

    _logger.info('prereduction ended')

    return result


class FiberFlatRecipe(MegaraBaseRecipe):
    '''Process FIBER_FLAT images and create MASTER_FIBER_FLAT.'''

    # Requirements
    master_bias = MasterBiasRequirement()
    obresult = ObservationResultRequirement()
    # Products
    fiberflat_frame = Product(MasterFiberFlat)
    fiberflat_rss = Product(MasterFiberFlat)
    traces = Product(TraceMapType)

    def __init__(self):
        super(FiberFlatRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, rinput):
        return self.process_base1(rinput.obresult, rinput.master_bias)

    
    def process_base1(self, obresult, master_bias):
        _logger.info('starting fiber flat reduction')

        reduced = process_common(self, obresult, master_bias)
        
        cstart = 2000
        step = 2
    
        tracemap = self.trace(reduced[0].data, cstart, step)
        
        rss = self.extract(reduced[0].data, tracemap)
        
        rss_norm = rss
        
        _logger.info('fiber flat reduction ended')

        result = self.create_result(fiberflat_frame=reduced,
                                    fiberflat_rss=fits.PrimaryHDU(rss_norm),
                                    traces=tracemap)

        return result

    def trace(self, data, cstart, step):


        # fit_traces = domefun(data, cstart=2000, hs=20)

        cstart = cstart
        hs = 3
        step1 = step
        background1 = 150.0
        npred = 1
        maxdis1 = 2.0

        _logger.info('find peaks in column %i', cstart)

        central_peaks = init_traces(data, center=cstart, hs=hs,
                                background=background1, npred=npred)

        

        tracelist = []
        if data.dtype.byteorder != '=':
            _logger.debug('byteswapping image')
            image2 = data.byteswap().newbyteorder()
        else:
            image2 = data
            
        _logger.info('trace peaks')
        for trace in central_peaks.values():

            mm = tracing(image2, x=cstart, y=trace.trace_f[0], step=step1,
                         hs=hs, background=background1, maxdis=maxdis1)

            pfit = numpy.polyfit(mm[:,0], mm[:,1], deg=5)
            
            tracelist.append({'fibid': trace.fibid, 'boxid': trace.boxid,
                              'start':0, 'stop':4095,
                              'fitparms': pfit.tolist()})

            #plt.title('cython version, fiber %i' % trace.fibid)
            #plt.plot(mm[:,0], mm[:,1], 'r*')
            #xpix = np.arange(0, 2000, 1)
            #p = numpy.poly1d(pfit)
            #plt.plot(xpix, p(xpix), 'b')
            #plt.show()    

        result_traces = TraceMap(tracelist)
                      
        return result_traces

    def extract_region(self, data, border1, border2, pesos, xpos):
        
        extend = (border1.min(), border2.max())
        extend_pix = (wcs_to_pix(extend[0]), wcs_to_pix(extend[1])+1)
        region = slice(extend_pix[0],extend_pix[1])

        for x, a,b in zip(xpos, border1, border2):
            fill_other(pesos[:,x], a, b)

            final2d = data[region,:] * pesos[region,:]

        pesos[region,:] = 0.0
        final = final2d.sum(axis=0)
        return final


    def extract(self, data, tracemap):

        pols = [np.poly1d(trace['fitparms']) for trace in tracemap.traces]

        rss = np.empty((len(pols), data.shape[1]))

        xpos = np.arange(data.shape[1])

        pesos = np.zeros_like(data)


        # extract the first trace
        p1 = pols[0]
        p2 = pols[1]
        pix_2 = p2(xpos)
        pix_1 = p1(xpos)
        pix_12 = 0.5 * (pix_2 + pix_1)
        # Use the half distance in the first trace
        pix_01 = 1.5 * pix_1 - 0.5 * pix_2
        # Extract on the other side
        rss[0] = self.extract_region(data, pix_01, pix_12, pesos, xpos)
        
        fibid = 2
        for p2 in pols[2:-1]:

            print ('fiber', fibid)            
            pix_1, pix_01 = pix_2, pix_12 
            pix_2 = p2(xpos)

            pix_12 = 0.5 * (pix_2 + pix_1)
            
            rss[fibid - 1] = self.extract_region(data, pix_01, pix_12, pesos, xpos)
            
            fibid += 1

        return rss

class TwilightFiberFlatRecipe(MegaraBaseRecipe):

    master_bias = MasterBiasRequirement()
    obresult = ObservationResultRequirement()

    fiberflat_frame = Product(MasterFiberFlat)
    fiberflat_rss = Product(MasterFiberFlat)
    traces = Product(ArrayType)

    def __init__(self):
        super(TwilightFiberFlatRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, rinput):
        pass


class TraceMapRecipe(MegaraBaseRecipe):

    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    fiberflat_frame = Product(MasterFiberFlat)
    traces = Product(TraceMapType)

    def __init__(self):
        super(TraceMapRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, rinput):

        result = self.process_base(rinput.obresult, rinput.master_bias)

        data = result[0].data

        # fit_traces = domefun(data, cstart=2000, hs=20)

        cstart = 2000
        hs = 3
        step1 = 2
        background1 = 150.0
        npred = 1
        maxdis1 = 2.0

        _logger.info('find peaks in column %i', cstart)

        central_peaks = init_traces(data, center=cstart, hs=hs,
                                background=background1, npred=npred)

        

        tracelist = []
        if data.dtype.byteorder != '=':
            _logger.debug('byteswapping image')
            image2 = data.byteswap().newbyteorder()
        else:
            image2 = data
            
        _logger.info('trace peaks')
        for trace in central_peaks.values():

            mm = tracing(image2, x=cstart, y=trace.trace_f[0], step=step1,
                         hs=hs, background=background1, maxdis=maxdis1)

            pfit = numpy.polyfit(mm[:,0], mm[:,1], deg=5)
            
            tracelist.append({'fibid': trace.fibid, 'boxid': trace.boxid,
                              'start':0, 'stop':4095,
                              'fitparms': pfit.tolist()})

            #plt.title('cython version, fiber %i' % trace.fibid)
            #plt.plot(mm[:,0], mm[:,1], 'r*')
            #xpix = np.arange(0, 2000, 1)
            #p = numpy.poly1d(pfit)
            #plt.plot(xpix, p(xpix), 'b')
            #plt.show()    

        result_traces = TraceMap(tracelist)
                      
        return self.create_result(fiberflat_frame=result,
                                  traces=result_traces)

    def process_base(self, obresult, master_bias):
        reduced = process_common(self, obresult, master_bias)
        return reduced
