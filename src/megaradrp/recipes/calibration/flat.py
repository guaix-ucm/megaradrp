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
from megaradrp.products import TraceMapType
from megaradrp.requirements import MasterBiasRequirement

from .traces import domefun


_logger = logging.getLogger('numina.recipes.megara')


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
    traces = Product(ArrayType)

    def __init__(self):
        super(FiberFlatRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, rinput):
        return self.process_base(rinput.obresult, rinput.master_bias)

    def process_base(self, obresult, master_bias):
        _logger.info('starting fiber flat reduction')

        reduced = process_common(self, obresult, master_bias)
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

        fit_traces = domefun(data, cstart=2000, hs=20)

        return self.create_result(fiberflat_frame=result,
                                  traces=fit_traces)

    def process_base(self, obresult, master_bias):
        reduced = process_common(self, obresult, master_bias)
        return reduced
