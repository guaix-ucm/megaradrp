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
from megaradrp.processing import OverscanCorrector, TrimImage
# from numina.logger import log_to_history

from megaradrp.products import MasterFiberFlat
from megaradrp.products import TraceMap
from megaradrp.requirements import MasterBiasRequirement

from megaradrp.core import apextract_tracemap

_logger = logging.getLogger('numina.recipes.megara')


class ArcCalibrationRecipe(MegaraBaseRecipe):
    '''Process ARC images and create WL_CALIBRATION.'''

    # Requirements
    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    master_trace = Requirement(TraceMap, 'Trace information of the Apertures')
    line_catalog = Requirement(ArrayType, 'Catalog of lines')
    # Products
    arc_image = Product(DataFrameType)
    arc_rss = Product(DataFrameType)

    def __init__(self):
        super(ArcCalibrationRecipe, self).__init__(
            version="0.1.0"
        )

    def run(self, rinput):
        # Basic processing
        reduced = self.process_common(rinput.obresult, rinput.master_bias)

        _logger.info('extract fibers')
        rss = apextract_tracemap(reduced[0].data, rinput.tracemap)

        # WL calibration goes here
        return self.create_result(arc_image=reduced, arc_rss=rss)

    def process_common(self, obresult, master_bias):
        _logger.info('starting prereduction')

        o_c = OverscanCorrector()
        t_i = TrimImage()

        with master_bias.open() as hdul:
            mbias = hdul[0].data.copy()
            b_c = BiasCorrector(mbias)

        basicflow = SerialFlow([o_c, t_i, b_c])

        cdata = []

        try:
            for frame in obresult.images:
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
