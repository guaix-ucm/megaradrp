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

"""Arc Calibration Recipe for Megara"""

from __future__ import division, print_function

import logging

import numpy
from astropy.io import fits
from scipy.interpolate import interp1d

from numina.core import Requirement, Product, Parameter
from numina.core import DataFrameType
from numina.core.products import ArrayType
from numina.core.requirements import ObservationResultRequirement
from numina.core.products import LinesCatalog
from numina.array.wavecal.arccalibration import arccalibration_direct
from numina.array.wavecal.arccalibration import fit_solution
from numina.array.wavecal.arccalibration import gen_triplets_master
from numina.array.wavecal.statsummary import sigmaG
# FIXME: still using findpeaks1D functions
from numina.array.peaks.findpeaks1D import findPeaks_spectrum
from numina.array.peaks.findpeaks1D import refinePeaks_spectrum

from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.products import TraceMap
from megaradrp.requirements import MasterBiasRequirement, MasterBPMRequirement
from megaradrp.requirements import MasterDarkRequirement
from megaradrp.core.processing import apextract_tracemap

_logger = logging.getLogger('numina.recipes.megara')


class ArcCalibrationRecipe(MegaraBaseRecipe):
    """Process ARC images and create WL_CALIBRATION."""

    # Requirements
    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_bpm = MasterBPMRequirement()
    tracemap = Requirement(TraceMap, 'Trace information of the Apertures')
    lines_catalog = Requirement(LinesCatalog, 'Catalog of lines')
    polynomial_degree = Parameter(5, 'Polynomial degree of arc calibration')
    # Products
    arc_image = Product(DataFrameType)
    arc_rss = Product(DataFrameType)
    master_wlcalib = Product(ArrayType)

    def __init__(self):
        super(ArcCalibrationRecipe, self).__init__("0.1.0")

    def run(self, rinput):
        # Basic processing
        parameters = self.get_parameters(rinput)
        reduced = self.bias_process_common(rinput.obresult, parameters)

        _logger.info('extract fibers')
        _logger.info('extract fibers, %i', len(rinput.tracemap))
        rssdata = apextract_tracemap(reduced[0].data, rinput.tracemap)
        rsshdu = fits.PrimaryHDU(rssdata, header=reduced[0].header)
        rss = fits.HDUList([rsshdu])

        _logger.info('extracted %i fibers', rssdata.shape[0])

        # Skip any other inputs for the moment
        coeff_table = self.calibrate_wl(rssdata, rinput.lines_catalog,
                                        rinput.polynomial_degree)

        # WL calibration goes here
        return self.create_result(arc_image=reduced, arc_rss=rss,
                                  master_wlcalib=coeff_table)

    def calibrate_wl(self, rss, lines_catalog, poldeg, times_sigma=50.0):
        #
        # read master table (TBM) and generate auxiliary parameters (valid for
        # all the slits) for the wavelength calibration
        wv_master = lines_catalog[:, 0]
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
            threshold = numpy.median(row) + times_sigma * sigmaG(row)
            ipeaks_int = findPeaks_spectrum(row, nwinwidth, threshold)
            # refine peaks fitting an appropriate function (providing float
            # numbers)
            ipeaks_float = refinePeaks_spectrum(row, ipeaks_int, nwinwidth)

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
                                                 poly_degree_wfit=poldeg,
                                                 times_sigma_polfilt=10.0,
                                                 times_sigma_inclusion=5.0)

                _logger.info('Solution for row %d completed', idx)
                _logger.info('Fitting solution for row %d', idx)
                numpy_array_with_coeff, crval1_approx, cdelt1_approx = \
                    fit_solution(wv_master,
                                 xpeaks_refined,
                                 solution,
                                 poly_degree_wfit=poldeg,
                                 weighted=False)

                _logger.info('approximate crval1, cdelt1: %f %f',
                             crval1_approx,
                             cdelt1_approx)

                _logger.info('fitted coefficients %s', numpy_array_with_coeff)

                coeff_table[idx] = numpy_array_with_coeff

            except TypeError as error:
                _logger.error("%s", error)

        return coeff_table
