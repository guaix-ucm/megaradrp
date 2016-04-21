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

"""Fiber flat calibration Recipe for Megara"""

from __future__ import division, print_function

import logging
import multiprocessing as mp
import numpy

from astropy.io import fits

from numina.core import Product, Requirement

from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.products import MasterFiberFlat
from megaradrp.products import WavelengthCalibration, MasterWeights
from megaradrp.requirements import MasterBiasRequirement, MasterBPMRequirement
from megaradrp.requirements import MasterDarkRequirement, MasterSlitFlatRequirement
from megaradrp.core.processing import apextract_tracemap
from numina.core.products import DataFrameType

_logger = logging.getLogger('numina.recipes.megara')


class FiberFlatRecipe(MegaraBaseRecipe):
    """Process FIBER_FLAT images and create MASTER_FIBER_FLAT."""

    # Requirements
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_bpm = MasterBPMRequirement()
    master_slitflat = MasterSlitFlatRequirement()
    wlcalib = Requirement(WavelengthCalibration, 'Wavelength calibration table')
    master_weights = Requirement(MasterWeights, 'Set of files')
    # Products
    fiberflat_frame = Product(DataFrameType)
    rss_fiberflat = Product(DataFrameType)
    master_fiberflat = Product(MasterFiberFlat)

    def __init__(self):
        super(FiberFlatRecipe, self).__init__(version="0.1.0")
        self.SIZE = 4096
        self.processes = mp.cpu_count()-2

    def run(self, rinput):
        # Basic processing
        parameters = self.get_parameters(rinput)

        _logger.info('process common')
        reduced = self.bias_process_common(rinput.obresult, parameters)

        # reduced.writeto('reduced_100.fits', clobber=True)
        _logger.info('Generada!')

        rss = numpy.array(self.load_files_from_directory(rinput.master_weights, img=reduced[0].data))

        _logger.info('Starting: resample_rss_flux')
        final, wcsdata = self.resample_rss_flux(rss.T, rinput.wlcalib)

        _logger.info('Compute mean and resampling again')

        mean = final[:,2000:2100].mean()
        aux = rss / mean

        # Add WCS spectral keywords
        hdu_f = fits.PrimaryHDU(aux.T, header=reduced[0].header)
        master_fiberflat = fits.HDUList([hdu_f])

        # fiberflat_frame = fits.PrimaryHDU(reduced, header=reduced[0].header)
        rss_fiberflat = fits.PrimaryHDU(final, header=reduced[0].header)

        result = self.create_result(master_fiberflat=master_fiberflat, fiberflat_frame=reduced, rss_fiberflat=rss_fiberflat)
        return result


    def resample_rss_flux(self, rss_old, wcalib):
        """Resample conserving the flux."""
        import math
        from numpy.polynomial.polynomial import polyval
        from numina.array.interpolation import SteffenInterpolator

        nfibers = rss_old.shape[0]
        nsamples = rss_old.shape[1]
        z = [0, nsamples-1]
        res = polyval(z, wcalib.T)
        all_delt = (res[:,1] - res[:,0]) / nsamples

        delts = all_delt.min()

        # first pixel is
        wl_min = res[:,0].min()
        # last pixel is
        wl_max = res[:,1].max()

        npix = int(math.ceil((wl_max - wl_min) / delts))
        new_x = numpy.arange(npix)
        new_wl = wl_min + delts * new_x

        old_x_borders = numpy.arange(-0.5, nsamples)
        old_wl_borders = polyval(old_x_borders, wcalib.T)

        new_borders = self.map_borders(new_wl)

        accum_flux = numpy.empty((nfibers, nsamples+1))
        accum_flux[:, 1:] = numpy.cumsum(rss_old, axis=1)
        accum_flux[:, 0] = 0.0
        rss_resampled = numpy.zeros((nfibers, npix))
        for idx in range(nfibers):
            # We need a monotonic interpolator
            # linear would work, we use a cubic interpolator
            interpolator = SteffenInterpolator(old_wl_borders[idx], accum_flux[idx], extrapolate='border')
            fl_borders = interpolator(new_borders)
            rss_resampled[idx] = fl_borders[1:]- fl_borders[:-1]
        return rss_resampled, (wl_min, wl_max, delts)

    def map_borders(self, wls):
        """Compute borders of pixels for interpolation.

        The border of the pixel is assumed to be midway of the wls
        """
        midpt_wl = 0.5 * (wls[1:] + wls[:-1])
        all_borders = numpy.zeros((wls.shape[0]+1,))
        all_borders[1:-1] = midpt_wl
        all_borders[0] = 2 * wls[0] - midpt_wl[0]
        all_borders[-1] = 2 * wls[-1] - midpt_wl[-1]
        return all_borders

    def add_wcs(self, hdr, wlr0, delt):
        hdr['CRPIX1'] = 1
        hdr['CRVAL1'] = wlr0
        hdr['CDELT1'] = delt
        hdr['CTYPE1'] = 'WAVELENGTH'
        hdr['CRPIX2'] = 1
        hdr['CRVAL2'] = 1
        hdr['CDELT2'] = 1
        hdr['CTYPE2'] = 'PIXEL'
        return hdr

    def decompress(self, tar_file):
        '''
        :param tar_name: <str> name of the tar file
        :return: None
        '''
        import tarfile

        name = tar_file.fileobj.name.split('.tar')[0]

        aux = tar_file.extractall(name+'/')
        return name

    def load_files_from_directory(self, tar_file, img):
        '''
        :param img: <numpy.array> reduced image
        :param tar_file: <pointer> tar file to be extracted
        :return: list of extracted weights
        '''

        path = self.decompress(tar_file)

        _logger.info('decompress done')
        _logger.info('Starting: _load_files_paralell')
        pool = mp.Pool(processes=self.processes)
        results = [pool.apply_async(_load_files_paralell,
                                    args=(ite, path)) for ite in range(self.SIZE)]
        results = [p.get() for p in results]
        # return results

        _logger.info('Starting: extract_w_paralell')

        pool2 = mp.Pool(processes=self.processes)
        extracted_w = [pool2.apply_async(extract_w_paralell,
                                    args=(img[:,ite], results[ite])) for ite in range(self.SIZE)]
        extracted_w = [p.get() for p in extracted_w]

        _logger.info('extracted')

        return extracted_w

    def run_tracemap(self, rinput):
        # Basic processing
        parameters = self.get_parameters(rinput)

        reduced = self.bias_process_common(rinput.obresult, parameters)

        _logger.info('extract fibers')
        rssdata = apextract_tracemap(reduced[0].data, rinput.tracemap)
        # FIXME: we are ignoring here all the possible bad pixels
        # and WL distortion when doing the normalization
        # rssdata /= rssdata.mean() #Originally uncomment
        rsshdu = fits.PrimaryHDU(rssdata, header=reduced[0].header)
        rss = fits.HDUList([rsshdu])

        _logger.info('extraction completed')
        _logger.info('fiber flat reduction ended')

        result = self.create_result(fiberflat_frame=reduced, master_fiberflat=rss)
        return result

def extract_w_paralell(img, mlist):

    '''
    :param img: <fits>
    :param mlist: <list> one element of the csr_matrix
    :return: <ndarray> result of lsqr
    '''
    from scipy.sparse.linalg import lsqr
    x = lsqr(mlist, img)
    return x[0]

def _load_files_paralell(col, path):
    '''
    :param col: <str,int> name of the fits file. It is a counter
    :param path: <str> path where *.npz are
    :return: csr_matrix
    '''
    from scipy.sparse import csr_matrix

    filename = '%s/%s.npz' % (path, col)
    loader = numpy.load(filename)
    return csr_matrix(
        (loader['data'], loader['indices'], loader['indptr']),
        shape=loader['shape'])