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

from numina.core import Requirement, Product, Parameter
from numina.core import DataFrameType
from numina.core.products import ArrayType
from numina.core.requirements import ObservationResultRequirement
from numina.array.combine import median as c_median
from numina.flow import SerialFlow
from numina.flow.processing import BiasCorrector

# For WL calibration
# FIXME: remove this later
from numina.core.products import LinesCatalog
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from numina.array.wavecal.arccalibration import arccalibration_direct
from numina.array.wavecal.arccalibration import fit_solution
from numina.array.wavecal.arccalibration import gen_triplets_master
from numina.array.wavecal.statsummary import sigmaG
from numina.array.peaks.findpeaks1D import findPeaks_spectrum
from numina.array.peaks.findpeaks1D import refinePeaks_spectrum

from megaradrp.core import MegaraBaseRecipe
from megaradrp.processing import OverscanCorrector, TrimImage

from megaradrp.products import TraceMap
from megaradrp.requirements import MasterBiasRequirement

from megaradrp.core import apextract_tracemap

_logger = logging.getLogger('numina.recipes.megara')


class ArcCalibrationRecipe(MegaraBaseRecipe):
    '''Process ARC images and create WL_CALIBRATION.'''

    # Requirements
    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    tracemap = Requirement(TraceMap, 'Trace information of the Apertures')
    lines_catalog = Requirement(LinesCatalog, 'Catalog of lines')
    polynomial_degree = Parameter(2, 'Polynomial degree of arc calibration')
    # Products
    arc_image = Product(DataFrameType)
    arc_rss = Product(DataFrameType)
    wlcalib = Product(ArrayType)

    def __init__(self):
        super(ArcCalibrationRecipe, self).__init__(
            version="0.1.0"
        )

    def run(self, rinput):
        # Basic processing
        reduced = self.process_common(rinput.obresult, rinput.master_bias)

        _logger.info('extract fibers')
        rssdata = apextract_tracemap(reduced[0].data, rinput.tracemap)
        rsshdu = fits.PrimaryHDU(rssdata, header=reduced[0].header)
        rss = fits.HDUList([rsshdu])

        # Skip any other inputs for the moment
        coeff_table = self.calibrate_wl(rssdata, rinput.lines_catalog,
                                        rinput.polynomial_degree)

        # WL calibration goes here
        return self.create_result(arc_image=reduced, arc_rss=rss,
                                  wlcalib=coeff_table)

    def calibrate_wl(self, rss, lines_catalog, poldeg, times_sigma=50.0):
        # 
        # read master table (TBM) and generate auxiliary parameters (valid for
        # all the slits) for the wavelength calibration
        wv_master = lines_catalog[:,0]
        ntriplets_master, ratios_master_sorted, triplets_master_sorted_list = \
                  gen_triplets_master(wv_master)
        # FIXME: this depends on the spectral and dispersion axes
        nspec = rss.shape[0]
        coeff_table = numpy.zeros((nspec, poldeg + 1))
        # Loop over rows in RSS
        nwinwidth = 5
        for idx, row in enumerate(rss):
            _logger.info('Starting row %d', idx)
            # find peaks (initial search providing integer numbers)
            threshold = numpy.median(row)+times_sigma*sigmaG(row)
            ipeaks_int = findPeaks_spectrum(row, nwinwidth=nwinwidth, 
                                                 data_threshold=threshold)
            # refine peaks fitting an appropriate function (providing float 
            # numbers)
            ipeaks_float = refinePeaks_spectrum(row, ipeaks_int, nwinwidth, 
                                                method=2)

   
            # define interpolation function and interpolate the refined peak 
            # location, passing from index number (within the row array) 
            # to channel number (note that this step takes care of the fact 
            # that the extracted spectrum may correspond to a subregion in the 
            # spectral direction)

            # FIXME: xchannel ???
            # This comes from Nico's code, so probably pixels
            # will start in 1
            naxis1 = row.shape[0]
            xchannel = numpy.arange(1, naxis1 + 1)

            finterp_channel = interp1d(range(xchannel.size), xchannel, 
                                       kind='linear')
            xpeaks_refined = finterp_channel(ipeaks_float)
    
            if False: # TBR (to be removed in the future)
                # plot extracted spectrum
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.set_title('row %d' % idx)
                ax.set_xlim([1,naxis1])
                ax.plot(xchannel,row,'k-')
                ax.plot(xchannel[ipeaks_int], row[ipeaks_int], 'ro')
                ax.plot([1,naxis1],[threshold, threshold], linestyle='dashed')
                ylim = ax.get_ylim()
                for xdum in zip(xpeaks_refined,xpeaks_refined):
                    ax.plot(xdum, ylim, linestyle='dotted', color='magenta')
                plt.show()
                #plt.show(block=False)
                #input('press <RETURN> to continue...')

           # wavelength calibration
            try:
                solution = arccalibration_direct(wv_master,
                                                 ntriplets_master,
                                                 ratios_master_sorted,
                                                 triplets_master_sorted_list,
                                                 xpeaks_refined,
                                                 naxis1,
                                                 wv_ini_search=3500,
                                                 wv_end_search=4500,
                                                 error_xpos_arc=2.0,
                                                 times_sigma_r=3.0,
                                                 frac_triplets_for_sum=0.50,
                                                 times_sigma_theil_sen=10.0,
                                                 poly_degree_wfit=2,
                                                 times_sigma_polfilt=10.0,
                                                 times_sigma_inclusion=5.0)

                _logger.info('Solution for row %d completed', idx)
                _logger.info('Fitting solution for row %d', idx)
                numpy_array_with_coeff, crval1_approx, cdelt1_approx = \
                  fit_solution(wv_master,
                               xpeaks_refined,
                               solution,
                               poly_degree_wfit=2,
                               weighted=False)
                
                _logger.info('approximate crval1, cdelt1: %f %f',
                             crval1_approx,cdelt1_approx)
                _logger.info('fitted coefficients %s',numpy_array_with_coeff)
                coeff_table[idx] = numpy_array_with_coeff
            except TypeError as error:
                _logger.error("%s", error)

        return coeff_table


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
