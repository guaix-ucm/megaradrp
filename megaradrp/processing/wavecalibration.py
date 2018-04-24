#
# Copyright 2016-2018 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Corrector for wavecalibration"""

import logging
import datetime

import numpy
from numpy.polynomial.polynomial import polyval
from astropy.io import fits
import numina.array.utils as u
from numina.processing import Corrector
from numina.array.interpolation import SteffenInterpolator

from megaradrp.instrument import WLCALIB_PARAMS

_logger = logging.getLogger(__name__)


class WavelengthCalibrator(Corrector):
    """A Node that applies wavelength calibration."""

    def __init__(self, solutionwl, datamodel=None, dtype='float32'):

        super(WavelengthCalibrator, self).__init__(
            datamodel=datamodel,
            calibid=solutionwl.calibid,
            dtype=dtype)

        self.solutionwl = solutionwl

    def run(self, rss):
        imgid = self.get_imgid(rss)
        _logger.debug('wavelength calibration in image %s', imgid)
        _logger.debug('with wavecalib %s', self.calibid)
        _logger.debug('offsets are %s', self.solutionwl.global_offset.coef)

        current_vph = rss[0].header['VPH']
        current_insmode = rss[0].header['INSMODE']

        _logger.debug('Current INSMODE is %s, VPH is %s', current_insmode, current_vph)
        if current_insmode in WLCALIB_PARAMS and current_vph in WLCALIB_PARAMS[current_insmode]:
            wvpar_dict = WLCALIB_PARAMS[current_insmode][current_vph]
            _logger.info('precomputed wl parameters are %s', wvpar_dict)
        else:
            msg = 'insmode {} grism {} is not defined in megaradrp.instrument.WLCALIB_PARAMS'.format(
                current_insmode,
                current_vph
            )
            raise ValueError(msg)

        _logger.debug('Resample RSS')
        final, wcsdata, values = self.resample_rss_flux(rss[0].data, wvpar_dict)

        _logger.debug('Update headers')
        rss_wl = fits.PrimaryHDU(
            data=final.astype(self.dtype),
            header=rss[0].header
        )

        hdr = rss_wl.header
        self.add_wcs(rss_wl.header, wvpar_dict['crval'], wvpar_dict['cdelt'],
                     wvpar_dict['crpix'])

        hdr['NUM-WAV'] = self.calibid
        hdr['history'] = 'Wavelength calibration with {}'.format(self.calibid)
        hdr['history'] = 'Aperture extraction offsets are {}'.format(
            self.solutionwl.global_offset.coef.tolist())
        hdr['history'] = 'Wavelength calibration time {}'.format(datetime.datetime.utcnow().isoformat())

        # Update other HDUs if needed
        rss[0] = rss_wl

        map_data = numpy.zeros(rss[0].shape, dtype='int16')

        fibers_ext = rss['FIBERS']
        fibers_ext_headers = fibers_ext.header

        for fibid, (s1, s2) in values:
            idx = fibid - 1
            map_data[idx, s1:s2+1] = 1
            # Update Fibers
            key = "FIB%03dW1" % fibid
            fibers_ext_headers[key] =  (s1 + 1, "Start of spectral coverage")
            key = "FIB%03dW2" % fibid
            fibers_ext_headers[key] =  (s2 + 1, "End of spectral coverage")

        for fibid in self.solutionwl.error_fitting:
            # Update Fibers
            key = "FIB%03d_V" % fibid
            fibers_ext_headers[key] =  False

        for fibid in self.solutionwl.missing_fibers:
            # Update Fibers
            key = "FIB%03d_V" % fibid
            fibers_ext_headers[key] =  False

        rss_map = fits.ImageHDU(data=map_data, name='WLMAP')

        rss.append(rss_map)
        return rss

    def add_wcs(self, hdr, wlr0, delt, crpix=0.0):
        c_crpix = 'Pixel coordinate of reference point'
        c_cunit = 'Units of coordinate increment and value'
        unit = 'Angstrom'
        c_crval = '[{}] Coordinate value at reference point'.format(unit)
        c_cdelt = '[{}] Coordinate increment at reference point'.format(unit)
        hdr['CRPIX1'] = (crpix, c_crpix)
        hdr['CRVAL1'] = wlr0, c_crval
        hdr['CDELT1'] = delt, c_cdelt
        hdr['CUNIT1'] = unit, c_cunit
        hdr['CTYPE1'] = 'WAVELENGTH'
        hdr['CRPIX2'] = (0.0, c_crpix)
        hdr['CRVAL2'] = 0.0
        hdr['CDELT2'] = 1.0
        hdr['CTYPE2'] = ''
        return hdr

    def resample_rss_flux(self, rss_old, wvpar_dict):

        nfibers = rss_old.shape[0]
        nsamples = rss_old.shape[1]

        npix = wvpar_dict['npix']
        delts = wvpar_dict['cdelt']
        wl_min = wvpar_dict['crval']
        crpix = wvpar_dict['crpix']

        wl_max = wl_min + (npix - crpix) * delts

        new_x = numpy.arange(npix)
        new_wl = wl_min + delts * new_x

        old_x_borders = numpy.arange(-0.5, nsamples)
        old_x_borders += crpix  # following FITS criterium
        # old_wl_borders = polyval(old_x_borders, wcalib.T)

        new_borders = self.map_borders(new_wl)

        accum_flux = numpy.empty((nfibers, nsamples + 1))
        accum_flux[:, 1:] = numpy.cumsum(rss_old, axis=1)
        accum_flux[:, 0] = 0.0
        rss_resampled = numpy.zeros((nfibers, npix))
        values = []

        for fibsol in self.solutionwl.contents:

            fibid = fibsol.fibid
            idx = fibid - 1

            # small correction defined in master_wlcalib_XXX_XX-X.json
            coeff = fibsol.solution.coeff
            offset_wavelength = self.solutionwl.global_offset(fibid)
            coeff[0] -= offset_wavelength

            old_wl_borders = polyval(old_x_borders, coeff)

            s1 = u.coor_to_pix_1d((old_wl_borders[0] - wl_min) / delts)
            s2 = u.coor_to_pix_1d((old_wl_borders[-1] - wl_min) / delts)

            s1 = max(0, min(s1, npix - 1))
            s2 = max(0, min(s2, npix - 1))

            # We need a monotonic interpolator
            # linear would work, we use a cubic interpolator
            interpolator = SteffenInterpolator(
                old_wl_borders,
                accum_flux[idx],
                extrapolate='border'
            )
            fl_borders = interpolator(new_borders)
            rss_resampled[idx] = fl_borders[1:] - fl_borders[:-1]
            values.append((fibid, (s1, s2)))

        return rss_resampled, (wl_min, wl_max, delts), values

    def map_borders(self, wls):
        """Compute borders of pixels for interpolation.

        The border of the pixel is assumed to be midway of the wls
        """
        midpt_wl = 0.5 * (wls[1:] + wls[:-1])
        all_borders = numpy.zeros((wls.shape[0] + 1,))
        all_borders[1:-1] = midpt_wl
        all_borders[0] = 2 * wls[0] - midpt_wl[0]
        all_borders[-1] = 2 * wls[-1] - midpt_wl[-1]
        return all_borders