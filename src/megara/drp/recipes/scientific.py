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

from numina import __version__
from numina.core import BaseRecipe, RecipeRequirements
from numina.core import Product, DataProductRequirement
from numina.core import define_requirements, define_result
from numina.core.requirements import ObservationResultRequirement, Requirement
from numina.array.combine import median as c_median
from numina.flow import SerialFlow
from numina.flow.processing import BiasCorrector

from megara.drp.core import OverscanCorrector, TrimImage
from megara.drp.core import ApertureExtractor, FiberFlatCorrector

#from numina.logger import log_to_history

from megara.drp.core import RecipeResult
from megara.drp.products import MasterBias, MasterFiberFlat
from megara.drp.products import MasterSensitivity,  TraceMapType
_logger = logging.getLogger('numina.recipes.megara')

class FiberMOSRecipeRequirements(RecipeRequirements):
    obresult = ObservationResultRequirement()
    master_bias = DataProductRequirement(MasterBias, 'Master bias calibration')
    master_fiber_flat = DataProductRequirement(MasterFiberFlat, 'Master fiber flat calibration')
    traces = Requirement(TraceMapType, 'Trace information of the Apertures')
    sensitivity = DataProductRequirement(MasterSensitivity, 'Sensitivity', optional=True)

class FiberMOSRecipeResult(RecipeResult):
    final = Product(MasterFiberFlat)
    target = Product(MasterFiberFlat)
    sky = Product(MasterFiberFlat)

@define_requirements(FiberMOSRecipeRequirements)
@define_result(FiberMOSRecipeResult)
class FiberMOSRecipe(BaseRecipe):
    '''Process MOS images.'''

    def __init__(self):
        super(FiberMOSRecipe, self).__init__(
                        author="Sergio Pascual <sergiopr@fis.ucm.es>",
                        version="0.1.0"
                )

    def run(self, rinput):
        _logger.info('starting fiber MOS reduction')
        
        
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
        s_data = []
        
        try:
            for frame in rinput.obresult.frames:
                hdulist = frame.open()
                hdulist = basicflow(hdulist)
                p_type = hdulist[0].header.get('OBSTYPE')
                if p_type == 'SKY':
                    s_data.append(hdulist)
                else:
                    t_data.append(hdulist)
 
            _logger.info('stacking %d sky images using median', len(s_data))
            
            data_s = c_median([d[0].data for d in s_data], dtype='float32')
            template_header = s_data[0][0].header
            hdu_s = fits.PrimaryHDU(data_s[0], header=template_header)
            
            data_t = c_median([d[0].data for d in t_data], dtype='float32')
            template_header = t_data[0][0].header
            hdu_t = fits.PrimaryHDU(data_t[0], header=template_header)
        finally:
            for hdulist in t_data:
                hdulist.close()
            for hdulist in s_data:
                hdulist.close()
      
        wlr = (3673.12731884058, 4417.497427536232)
        size = hdu_t.data.shape[1]    
        delt = (wlr[1] - wlr[0]) / (size-1)
        
        def add_wcs(hdr, numtyp):
            hdr['CRPIX1'] = 1
            hdr['CRVAL1'] = wlr[0]
            hdr['CDELT1'] = delt
            hdr['CTYPE1'] = 'WAVELENGTH'
            hdr['CRPIX2'] = 1
            hdr['CRVAL2'] = 1
            hdr['CDELT2'] = 1
            hdr['CTYPE2'] = 'PIXEL'
            hdr['NUMXVER'] = (__version__, 'Numina package version')
            hdr['NUMRNAM'] = (self.__class__.__name__, 'Numina recipe name')
            hdr['NUMRVER'] = (self.__version__, 'Numina recipe version')
            hdr['CCDMEAN'] = data_s[0].mean()
            hdr['NUMTYP'] = (numtyp, 'Data product type')
            return hdr
      
        add_wcs(hdu_s.header, 'SCIENCE_SKY')
        
        add_wcs(hdu_t.header, 'SCIENCE_TARGET')
            
        _logger.info('subtract SKY RSS from target RSS')
        final = data_t[0] - data_s[0]
        hdu_f = fits.PrimaryHDU(final, header=template_header)
                
        add_wcs(hdu_f.header, 'SCIENCE_FINAL')
                
        if rinput.sensitivity:
            _logger.info('apply sensitivity')
            with rinput.sensitivity.open() as hdul:
                sens = hdul[0].data
                hdu_f.data *= sens
        else:
            _logger.info('sensitivity is not defined, ignoring')
        
        _logger.info('MOS reduction ended')

        result = FiberMOSRecipeResult(final=hdu_f, target=hdu_t, sky=hdu_s)
        return result

