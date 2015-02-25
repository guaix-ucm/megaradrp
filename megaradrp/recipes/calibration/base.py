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


from numina.core import Product, DataProductRequirement, Requirement
from numina.core.products import ArrayType
from numina.core.requirements import ObservationResultRequirement
from numina.core import RecipeError
from numina.array.combine import median as c_median
from numina.flow import SerialFlow
from numina.flow.processing import BiasCorrector

from megaradrp.core import MegaraBaseRecipe
from megaradrp.core import OverscanCorrector, TrimImage
from megaradrp.core import ApertureExtractor, FiberFlatCorrector
from megaradrp.core import peakdet
# from numina.logger import log_to_history
from megaradrp.requirements import MasterBiasRequirement
from megaradrp.requirements import MasterFiberFlatRequirement
from megaradrp.products import MasterBias, MasterDark, MasterFiberFlat
from megaradrp.products import TraceMapType, MasterSensitivity


_logger = logging.getLogger('numina.recipes.megara')


class BiasRecipe(MegaraBaseRecipe):
    '''Process BIAS images and create MASTER_BIAS.'''

    obresult = ObservationResultRequirement()

    biasframe = Product(MasterBias)

    def __init__(self):
        super(BiasRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, rinput):
        return self.process(rinput.obresult)

    def process(self, obresult):
        _logger.info('starting bias reduction')

        if not obresult.frames:
            raise RecipeError('Frame list is empty')

        cdata = []

        o_c = OverscanCorrector()
        t_i = TrimImage()

        basicflow = SerialFlow([o_c, t_i])

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
        hdr = self.set_base_headers(hdr)
        hdr['IMGTYP'] = ('BIAS', 'Image type')
        hdr['NUMTYP'] = ('MASTER_BIAS', 'Data product type')
        hdr['CCDMEAN'] = data[0].mean()

        varhdu = fits.ImageHDU(data[1], name='VARIANCE')
        num = fits.ImageHDU(data[2], name='MAP')
        hdulist = fits.HDUList([hdu, varhdu, num])
        _logger.info('bias reduction ended')

        result = self.create_result(biasframe=hdu)
        return result


class DarkRecipe(MegaraBaseRecipe):

    '''Process DARK images and provide MASTER_DARK. '''

    master_bias = MasterBiasRequirement()

    darkframe = Product(MasterDark)

    def __init__(self):
        super(DarkRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    # FIXME find a better way of doing this automatically
    # @log_to_history(_logger)
    def run(self, rinput):

        _logger.info('starting dark reduction')

        _logger.info('dark reduction ended')

        result = self.create_result(darkframe=None)
        return result


class PseudoFluxCalibrationRecipe(MegaraBaseRecipe):

    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    master_fiber_flat = MasterFiberFlatRequirement()
    traces = Requirement(TraceMapType, 'Trace information of the Apertures')
    reference_spectrum = DataProductRequirement(
        MasterFiberFlat, 'Reference spectrum')

    calibration = Product(MasterSensitivity)
    calibration_rss = Product(MasterSensitivity)

    def __init__(self):
        super(PseudoFluxCalibrationRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, rinput):
        _logger.info('starting pseudo flux calibration')

        o_c = OverscanCorrector()
        t_i = TrimImage()

        with rinput.master_bias.open() as hdul:
            mbias = hdul[0].data.copy()
            b_c = BiasCorrector(mbias)

        a_e = ApertureExtractor(rinput.traces)

        with rinput.master_fiber_flat.open() as hdul:
            f_f_c = FiberFlatCorrector(hdul)

        basicflow = SerialFlow([o_c, t_i, b_c, a_e, f_f_c])

        t_data = []

        try:
            for frame in rinput.obresult.frames:
                hdulist = frame.open()
                hdulist = basicflow(hdulist)
                t_data.append(hdulist)

            data_t = c_median([d[0].data for d in t_data], dtype='float32')
            template_header = t_data[0][0].header
            hdu_t = fits.PrimaryHDU(data_t[0], header=template_header)
        finally:
            for hdulist in t_data:
                hdulist.close()

        hdr = hdu_t.header
        hdr = self.set_base_headers(hdr)
        hdr['CCDMEAN'] = data_t[0].mean()
        hdr['NUMTYP'] = ('SCIENCE_TARGET', 'Data product type')

        # FIXME: hardcoded calibration
        # Polynomial that translates pixels to wl
        _logger.warning('using hardcoded LR-U spectral calibration')
        wlcal = [7.12175997e-10, -9.36387541e-06,
                 2.13624855e-01, 3.64665269e+03]
        plin = numpy.poly1d(wlcal)
        wl_n_r = plin(range(1, hdu_t.data.shape[1] + 1))  # Non-regular WL

        _logger.info('resampling reference spectrum')

        wlr = [3673.12731884058, 4417.497427536232]
        size = hdu_t.data.shape[1]
        delt = (wlr[1] - wlr[0]) / (size - 1)

        def add_wcs(hdr):
            hdr['CRPIX1'] = 1
            hdr['CRVAL1'] = wlr[0]
            hdr['CDELT1'] = delt
            hdr['CTYPE1'] = 'WAVELENGTH'
            hdr['CRPIX2'] = 1
            hdr['CRVAL2'] = 1
            hdr['CDELT2'] = 1
            hdr['CTYPE2'] = 'PIXEL'
            return hdr

        with rinput.reference_spectrum.open() as hdul:
            # Needs resampling
            data = hdul[0].data
            w_ref = wcs.WCS(hdul[0].header)
            # FIXME: Hardcoded values
            # because we do not have WL calibration
            pix = range(1, len(data) + 1)
            wl = w_ref.wcs_pix2world(pix, 1)
            # The 0 mean 0-based
            si = interp1d(wl, data)
            # Reference spectrum evaluated in the irregular WL grid
            final = si(wl_n_r)

        sens_data = final / hdu_t.data
        hdu_sens = fits.PrimaryHDU(sens_data, header=hdu_t.header)

        # Very simple wl calibration
        # add_wcs(hdu_sens.header)

        # add_wcs(hdu_t.header)

        _logger.info('pseudo flux calibration reduction ended')

        result = self.create_result(
            calibration=hdu_sens, calibration_rss=hdu_t)
        return result


class ArcRecipe(MegaraBaseRecipe):

    master_bias = MasterBiasRequirement()
    obresult = ObservationResultRequirement()

    fiberflat_frame = Product(MasterFiberFlat)
    fiberflat_rss = Product(MasterFiberFlat)
    traces = Product(ArrayType)

    def __init__(self):
        super(ArcRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, rinput):
        pass


class LCB_IFU_StdStarRecipe(MegaraBaseRecipe):

    master_bias = MasterBiasRequirement()
    obresult = ObservationResultRequirement()

    fiberflat_frame = Product(MasterFiberFlat)
    fiberflat_rss = Product(MasterFiberFlat)
    traces = Product(ArrayType)

    def __init__(self):
        super(LCB_IFU_StdStarRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, rinput):
        pass


class FiberMOS_StdStarRecipe(MegaraBaseRecipe):

    master_bias = MasterBiasRequirement()
    obresult = ObservationResultRequirement()

    fiberflat_frame = Product(MasterFiberFlat)
    fiberflat_rss = Product(MasterFiberFlat)
    traces = Product(ArrayType)

    def __init__(self):
        super(FiberMOS_StdStarRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, rinput):
        pass


class SensitivityFromStdStarRecipe(MegaraBaseRecipe):

    master_bias = MasterBiasRequirement()
    obresult = ObservationResultRequirement()

    fiberflat_frame = Product(MasterFiberFlat)
    fiberflat_rss = Product(MasterFiberFlat)
    traces = Product(ArrayType)

    def __init__(self):
        super(SensitivityFromStdStarRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, rinput):
        pass


class S_And_E_FromStdStarsRecipe(MegaraBaseRecipe):

    master_bias = MasterBiasRequirement()
    obresult = ObservationResultRequirement()

    fiberflat_frame = Product(MasterFiberFlat)
    fiberflat_rss = Product(MasterFiberFlat)
    traces = Product(ArrayType)

    def __init__(self):
        super(SensitivityFromStdStarRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, rinput):
        pass


class BadPixelsMaskRecipe(MegaraBaseRecipe):

    '''Process BIAS images and create MASTER_BIAS.'''

    obresult = ObservationResultRequirement()

    biasframe = Product(MasterBias)

    def __init__(self):
        super(BadPixelsMaskRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, rinput):
        pass


class LinearityTestRecipe(MegaraBaseRecipe):

    '''Process BIAS images and create MASTER_BIAS.'''

    obresult = ObservationResultRequirement()

    biasframe = Product(MasterBias)

    def __init__(self):
        super(LinearityTestRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, rinput):
        pass
