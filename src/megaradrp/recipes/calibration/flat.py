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

import logging

import numpy

from scipy.interpolate import interp1d

from astropy.io import fits
from astropy import wcs

from numina.core import BaseRecipeAutoQC
from numina.core import Product, DataProductRequirement, Requirement
from numina.core.products import ArrayType
from numina.core.requirements import ObservationResultRequirement
from numina.core import RecipeError
from numina.array.combine import median as c_median
from numina.flow import SerialFlow
from numina.flow.processing import BiasCorrector

from megaradrp.core import OverscanCorrector, TrimImage
from megaradrp.core import ApertureExtractor, FiberFlatCorrector
from megaradrp.core import peakdet
# from numina.logger import log_to_history

from megaradrp.products import MasterBias, MasterDark, MasterFiberFlat
from megaradrp.products import TraceMapType, MasterSensitivity

_logger = logging.getLogger('numina.recipes.megara')


class Trace(object):
    def __init__(self, start=0):
        self.start = start
        self.lost = None
        self.sample_c = []
        self.trace_c = []
        self.peak_c = []
        self.sample_f = []
        self.trace_f = []
        self.peak_f = []
        
    def predict_position(self, col):
        return self.trace_f[-1]


class FiberTrace(Trace):
    def __init__(self, fibid, start=0):
        super(FiberTrace, self).__init__(start=start)
        self.fibid = fibid


def delicate_centre(x, y):
    pos = numpy.polyfit(x, y, deg=2)
    tx = -pos[1] / (2 * pos[0])
    py = numpy.polyval(pos, tx)
    return tx, py, pos


class FiberFlatRecipe(BaseRecipeAutoQC):
    '''Process FIBER_FLAT images and create MASTER_FIBER_FLAT.'''

    # Requirements
    master_bias = DataProductRequirement(MasterBias, 'Master bias calibration')
    obresult = ObservationResultRequirement()
    # Products
    fiberflat_frame = Product(MasterFiberFlat)
    fiberflat_rss = Product(MasterFiberFlat)
    traces = Product(ArrayType)

    def __init__(self):
        super(FiberFlatRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, rinput):
        return self.process_base(rinput.obresult, rinput.master_bias)

    def process_common(self, obresult, master_bias):
        _logger.info('starting prep reduction')

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
        hdr = self.set_base_headers(hdr)
        hdr['CCDMEAN'] = data[0].mean()

        varhdu = fits.ImageHDU(data[1], name='VARIANCE')
        num = fits.ImageHDU(data[2], name='MAP')
        result = fits.HDUList([hdu, varhdu, num])

        _logger.info('prereduction ended')

        return result

    def process_base(self, obresult, master_bias):
        _logger.info('starting fiber flat reduction')

        reduced = self.process_common(obresult, master_bias)
        mm = reduced[0].data

        # Trace extract and normalize
        # Cut a region in the center

        cut = mm[:, 1980:2020]
        colcut = cut.sum(axis=1) / 40.0

        # Find peaks
        maxt, mint = peakdet(v=colcut, delta=0.3, back=5e3)
        _logger.info('found %d peaks', len(maxt))
        # Cut around the peak

        # Maximum half width of peaks
        maxw = 3

        borders = numpy.empty((maxt.shape[0], 3), dtype='int')
        borders[:, 1] = maxt[:, 0]
        borders[1:, 0] = mint[:-1, 0]
        borders[0, 0] = 0
        borders[:-1, 2] = mint[1:, 0]
        borders[-1, 2] = 1e4

        borders[:, 2] = numpy.minimum(borders[:, 2], borders[:, 1] + maxw)
        borders[:, 0] = numpy.maximum(borders[:, 0], borders[:, 1] - maxw)

        _logger.info('extract fibers')
        rss = numpy.empty((borders.shape[0], mm.shape[1]))
        for idx, r in enumerate(borders):
            l = int(r[0])
            r = int(r[2]) + 1
            sl = (slice(l, r), )
            m = mm[sl].sum(axis=0)
            rss[idx] = m

        # Normalize RSS
        _logger.info('normalize fibers')
        nidx = 350
        single_s = rss[nidx, :]

        # FIXME: hardcoded values
        x_fit_min = 50
        x_fit_max = 4050  # This is the last value, not included
        x_fit = range(x_fit_min, x_fit_max)

        # fitting a 3rd degree polynomial
        pol_coeff = numpy.polyfit(x_fit, single_s[x_fit_min:x_fit_max], deg=3)
        pol_fit = numpy.poly1d(pol_coeff)
        norm = pol_fit(range(single_s.shape[0]))
        rss_norm = rss / norm

        _logger.info('fiber flat reduction ended')

        result = self.create_result(fiberflat_frame=reduced,
                                    fiberflat_rss=fits.PrimaryHDU(rss_norm),
                                    traces=borders)

        return result

    def process_advanced(self, obresult, master_bias):
        _logger.info('starting fiber flat reduction')
        import matplotlib.pyplot as plt
        from .peakdetection import peak_detection_mean_window
        reduced = self.process_common(obresult, master_bias)
        image = reduced[0].data

        maxdis = 1.5
        xmin  = 1
        xmax = 4000
        ixmin= 100
        ixmax = 4000
        xplot = numpy.arange(xmin, xmax, 0.1)
        x = numpy.arange(xmin, xmax)
        data = image[xmin:xmax]
        #x = range(len(data))
        hs = 5

        # Number of steps needed to predict the position
        npredict = 0

        background = 0.0
        cstart = 2000
        cend = 50
        c = cstart
        # Trace extract and normalize
        # Cut a region in the center
        _logger.info('find groups of fibers')
        _logger.info('end find groups of fibers')
        cut_region = slice(c-hs, c+hs)
        cut = data[:,cut_region]
        colcut = cut.mean(axis=1)
        maxt = peak_detection_mean_window(colcut, x=x, k=3, xmin=ixmin, xmax=ixmax, background=background)
        npeaks = len(maxt)
        peakdist = numpy.diff(maxt[:,1])
        # number of peaks
        print 'npeaks', npeaks
        print 'distances between peaks (max)', peakdist.max(), '(min)', peakdist.min()
        #
        fiber_traces = {i: FiberTrace(i, start=c) for i in range(npeaks)}
        #
        plt.plot(x, colcut, 'r*-')
        plt.gca().set_title('npeaks =%i, c=%i, hs=%i' % (npeaks, c, hs))
        plt.scatter(maxt[:,1], 1.1*maxt[:,2], c='g')
        plt.savefig('fig1.png')
        # peaks
    
        for fibid, trace in fiber_traces.items():
            trace.sample_c.append(c)
            trace.trace_c.append(maxt[fibid,1])
            trace.peak_c.append(maxt[fibid,2])
            pixmax = int(maxt[fibid,0])
            # Take 2*2+1 pix
            tx, py, pos = delicate_centre(x[pixmax-2: pixmax+2+1], colcut[pixmax-2: pixmax+2+1])
            trace.sample_f.append(c)       
            trace.trace_f.append(tx)
            trace.peak_f.append(py)


        _logger.info('fiber flat reduction ended')

        result = self.create_result(fiberflat_frame=reduced,
                                    fiberflat_rss=None,
                                    traces=None)

        return result


class TwiligthFiberFlatRecipe(BaseRecipeAutoQC):

    master_bias = DataProductRequirement(MasterBias, 'Master bias calibration')
    obresult = ObservationResultRequirement()

    fiberflat_frame = Product(MasterFiberFlat)
    fiberflat_rss = Product(MasterFiberFlat)
    traces = Product(ArrayType)

    def __init__(self):
        super(TwiligthFiberFlatRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, rinput):
        pass


class TraceMapRecipe(BaseRecipeAutoQC):

    obresult = ObservationResultRequirement()
    biasframe = Product(MasterBias)

    def __init__(self):
        super(TraceMapRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, rinput):
        pass
