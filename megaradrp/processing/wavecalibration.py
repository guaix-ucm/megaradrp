#
# Copyright 2016-2019 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Corrector for wavelength calibration"""

import logging
import datetime
import warnings
import uuid

import numpy
from numpy.polynomial.polynomial import polyval
import astropy.wcs
import astropy.io.fits as fits

import numina.datamodel as dm
import numina.array.utils as utils
from numina.frame.utils import copy_img
from numina.processing import Corrector
from numina.array.interpolation import SteffenInterpolator

from megaradrp.instrument import WLCALIB_PARAMS

_logger = logging.getLogger(__name__)


class SimpleWcs1D(object):
    """Store parameters of a simple 1D WCS"""
    def __init__(self, crval=0, crpix=0, cdelt=1, ctype='', cunit='', **kwargs):
        self.crval = crval
        self.crpix = crpix
        self.cdelt = cdelt
        self.ctype = ctype
        self.cunit = cunit

    def create_internal_wcs_(self):
        # Not using units, only when writing in the header
        w = astropy.wcs.WCS(naxis=2, fix=False)
        w.wcs.crpix = [self.crpix, 0.0]
        w.wcs.cdelt = [self.cdelt, 1.0]
        w.wcs.crval = [self.crval, 0.0]
        w.wcs.ctype = ['AWAV', ' ']
        return w


class WavelengthCalibrator(Corrector):
    """A Node that applies wavelength calibration."""

    def __init__(self, solutionwl, datamodel=None, dtype='float32'):

        super(WavelengthCalibrator, self).__init__(
            datamodel=datamodel,
            calibid=solutionwl.calibid,
            dtype=dtype)

        self.solutionwl = solutionwl

    def run(self, rss):

        newrss = calibrate_wl_rss_megara(
            rss, self.solutionwl,
            dtype=self.dtype, span=2, inplace=True
        )
        return newrss


def calibrate_wl_rss_megara(rss, solutionwl, dtype='float32', span=0, inplace=False):
    """Apply wavelength calibration to a RSS

    Parameters
    ----------

    rss: astropy.io.fits.HDUList
        A Row stacked Spectra MEGARA image, not WL calibrated
    solutionwl: megaradrp.products.wavecalibration.WavelengthCalibration
        A wavelength calibration solution
    dtype: str
        Convertible to numpy.dtype
    span: int
        Remove `span` pixels at both sides of the resampled image
    inplace: bool
        The input RSS is modified inplace, if False, a copy is made

    Returns
    -------
    rss: astropy.io.fits.HDUList
        A Row stacked Spectra MEGARA image, WL calibrated
    """

    imgid = dm.get_imgid(rss)
    _logger.debug('wavelength calibration in image %s', imgid)
    _logger.debug('with wavecalib %s', solutionwl.calibid)
    _logger.debug('offsets are %s', solutionwl.global_offset.coef)

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

    targetwcs = SimpleWcs1D(**wvpar_dict)
    npix = wvpar_dict['npix']

    result = calibrate_wl_rss(
        rss, solutionwl, npix, targetwcs,
        dtype=dtype,
        span=span, inplace=inplace
    )
    return result


def calibrate_wl_rss(rss, solutionwl, npix, targetwcs, dtype='float32', span=0, inplace=False):
    """Apply wavelength calibration to a RSS

    Parameters
    ----------

    rss: astropy.io.fits.HDUList
        A Row stacked Spectra MEGARA image, not WL calibrated
    solutionwl: megaradrp.products.wavecalibration.WavelengthCalibration
        A wavelength calibration solution
    npix: int
        Number of channels of the calibrated RSS
    targetwcs: SimpleWcs1D
        Common WCS solution
    dtype: str
        Convertible to numpy.dtype
    span: int
        Remove `span` pixels at both sides of the resampled image
    inplace: bool
        The input RSS is modified inplace, if False, a copy is made

    Returns
    -------
    rss: astropy.io.fits.HDUList
        A Row stacked Spectra MEGARA image, WL calibrated
    """

    if not inplace:
        # This is a new HDUList
        rss = copy_img(rss)

    imgid = dm.get_imgid(rss)
    _logger.debug('wavelength calibration in image %s', imgid)
    _logger.debug('with wavecalib %s', solutionwl.calibid)
    _logger.debug('offsets are %s', solutionwl.global_offset.coef)

    # Target WCS
    re_wcs = targetwcs.create_internal_wcs_()

    _logger.debug('Resample RSS')
    final, limits = resample_rss_flux(
        rss[0].data, solutionwl, npix, re_wcs,
        span=span, fill=0
    )

    rss[0].data = final.astype(dtype)

    hdr = rss[0].header
    _logger.debug('Add WCS headers')
    rss_add_wcs(hdr, targetwcs.crval, targetwcs.cdelt, targetwcs.crpix)
    try:
        header_add_barycentric_correction(hdr, key='B')
    except KeyError as error:
        _logger.warning('Missing key %s, cannot add barycentric correction', error)
    _logger.debug('Add calibration headers')
    hdr['NUM-WAV'] = solutionwl.calibid
    hdr['history'] = 'Wavelength calibration with {}'.format(solutionwl.calibid)
    hdr['history'] = 'Aperture extraction offsets are {}'.format(
        solutionwl.global_offset.coef.tolist())
    hdr['history'] = 'Wavelength calibration time {}'.format(datetime.datetime.utcnow().isoformat())
    hdr['history'] = 'Resample span={}'.format(span)
    # Update UUID
    hdr['UUID'] = str(uuid.uuid1())

    # Update other HDUs if needed
    # dtype here can be int16 or uint8
    map_data = numpy.zeros_like(final, dtype='int16')

    fibers_ext = rss['FIBERS']
    fibers_ext_headers = fibers_ext.header
    # Add KEYWORDS
    # FIB%03dW1, FIB%03dW2
    for fibid, (lower, upper) in limits:
        idx = fibid - 1
        map_data[idx, lower:upper+1] = 1
        # Update Fibers
        key = "FIB{:03d}W1".format(fibid)
        fibers_ext_headers[key] =  (lower + 1, "Start of spectral coverage")
        key = "FIB{:03d}W2".format(fibid)
        fibers_ext_headers[key] =  (upper + 1, "End of spectral coverage")

    # Update KEYWORDS
    # "FIB%03d_V"
    for fibid in solutionwl.error_fitting:
        # Update Fibers
        key = "FIB%03d_V" % fibid
        fibers_ext_headers[key] =  False

    for fibid in solutionwl.missing_fibers:
        # Update Fibers
        key = "FIB%03d_V" % fibid
        fibers_ext_headers[key] =  False

    rss_map = fits.ImageHDU(data=map_data, name='WLMAP')

    rss.append(rss_map)
    return rss


def rss_add_wcs(hdr, crval, cdelt, crpix):
    """Add MEGARA 2D wavelength calibration headers"""
    c_crpix = 'Pixel coordinate of reference point'
    c_cunit = 'Units of coordinate increment and value'
    unit = 'Angstrom'
    c_crval = 'Coordinate value at reference point'
    c_cdelt = 'Coordinate increment at reference point'
    hdr['CRPIX1'] = crpix, c_crpix
    hdr['CRVAL1'] = crval, c_crval
    hdr['CDELT1'] = cdelt, c_cdelt
    hdr['CUNIT1'] = unit, c_cunit
    hdr['CTYPE1'] = 'AWAV', 'Air wavelength (linear)'

    hdr['CRPIX2'] = 0.0, c_crpix
    hdr['CRVAL2'] = 0.0, c_crval
    hdr['CDELT2'] = 1.0, c_cdelt
    hdr['CTYPE2'] = ''
    return hdr


def create_internal_wcs_(wlr0, delt, crpix):
    # Not using units, only when writing in the header
    w = astropy.wcs.WCS(naxis=2, fix=False)
    w.wcs.crpix = [crpix, 0.0]
    w.wcs.cdelt = [delt, 1.0]
    w.wcs.crval = [wlr0, 0.0]
    w.wcs.ctype = ['AWAV', ' ']

    return w


def header_add_barycentric_correction(hdr, key='B', out=None):
    """Add WCS keywords with barycentric correction

    Raises
    ------
    KeyError
        If a required keyword is missing
    TypeError
        If the header does not contain a spectral axis
    """
    from astropy.coordinates import SkyCoord, EarthLocation
    import astropy.time
    import astropy.constants as cons

    # Header must have DATE-OBS
    if 'DATE-OBS' not in hdr:
        raise KeyError("Keyword 'DATE-OBS' not found.")
    # Header must contain a primary WCS
    # Header must contain RADEG and DECDEG

    if 'OBSGEO-X' not in hdr:
        warnings.warn('OBSGEO- keywords not defined, using default values for GTC', RuntimeWarning)
        # Geocentric coordinates of GTC
        hdr['OBSGEO-X'] = 5327285.0921
        hdr['OBSGEO-Y'] = -1718777.1125
        hdr['OBSGEO-Z'] = 3051786.7327

    # Get main WCS
    wcs0 = astropy.wcs.WCS(hdr)
    if wcs0.wcs.spec == -1:
        # We don't have a spec axis
        raise TypeError('Header does not contain spectral axis')
    gtc = EarthLocation.from_geocentric(wcs0.wcs.obsgeo[0], wcs0.wcs.obsgeo[1], wcs0.wcs.obsgeo[2], unit='m')
    date_obs = astropy.time.Time(wcs0.wcs.dateobs, format='fits')
    # if frame='fk5', we need to pass the epoch and equinox
    sc = SkyCoord(ra=hdr['RADEG'], dec=hdr['DECDEG'], unit='deg')
    rv = sc.radial_velocity_correction(obstime=date_obs, location=gtc)
    factor = (1 + rv / cons.c).to('').value

    if out is None:
        out = hdr

    out['WCSNAME{}'.format(key)] = 'Barycentric correction'
    # out['CNAME1{}'.format(key)] = 'AxisV'
    out['CTYPE1{}'.format(key)] = hdr['CTYPE1']
    out['CRPIX1{}'.format(key)] = hdr['CRPIX1']
    out['CRVAL1{}'.format(key)] = hdr['CRVAL1'] * factor
    out['CDELT1{}'.format(key)] = hdr['CDELT1'] * factor
    out['CUNIT1{}'.format(key)] = hdr['CUNIT1']

    for keyword in ['CRPIX2', 'CRVAL2', 'CDELT2', 'CTYPE2']:
        try:
            out['{}{}'.format(keyword, key)] = hdr['{}'.format(keyword)]
        except KeyError:
            # Ignore non-existing key
            pass

    out['VELOSYS{}'.format(key)] = rv.to('m / s').value
    out['SPECSYS{}'.format(key)] = 'BARYCENT'
    out['SSYSOBS{}'.format(key)] = 'TOPOCENT'
    return out


def resample_rss_flux(arr, solutionwl, npix, finalwcs, span=0, fill=0):
    """Resample array according to a wavelength calibration solution

    Parameters
    ----------

    arr: numpy.ndarray
        A Row stacked Spectra MEGARA image, not WL calibrated
    solutionwl: megaradrp.products.wavecalibration.WavelengthCalibration
        A wavelength calibration solution
    npix: int
        Number of channels of the calibrated RSS
    finalwcs: astropy.wcs.WCS
        WCS solution of the final array
    span: int
        Remove `span` pixels at both sides of the resampled image
    fill: int
        Value used to fill the values removed by `span`

    Returns
    -------
    rss_resampled:  numpy.ndarray
        A resample array
    limits: a list of tuples
        Contains the fiberid and a pair with the first and last valid pixel (0-based)
    """

    nfibers = arr.shape[0]
    nsamples = arr.shape[1]

    # Use only the spectral axis
    subwcs = finalwcs.sub(['spectral'])

    # wl_max,  = subwcs.all_pix2world([npix], 1)

    # 0-based index
    new_x = numpy.arange(npix)
    new_wl,  = subwcs.all_pix2world(new_x, 0)

    # 0-based left borders
    old_x_borders_0 = numpy.arange(-0.5, nsamples)
    # 1-based left borders
    old_x_borders_1 = old_x_borders_0 + 1.0  # following FITS criterium

    # In AA
    new_wl_borders = pixel_borders(new_wl)

    accum_flux = numpy.empty((nfibers, nsamples + 1))
    accum_flux[:, 1:] = numpy.cumsum(arr, axis=1)
    accum_flux[:, 0] = 0.0
    rss_resampled = numpy.zeros((nfibers, npix))
    limits = []

    for fibsol in solutionwl.contents:

        fibid = fibsol.fibid
        idx = fibid - 1

        # small correction defined in master_wlcalib_XXX_XX-X.json
        coeff = fibsol.solution.coeff
        offset_wavelength = solutionwl.global_offset(fibid)
        coeff[0] -= offset_wavelength

        # Polynomial returns AA
        old_wl_borders = polyval(old_x_borders_1, coeff)
        # 0-based, AA
        ss_vals, = subwcs.all_world2pix(old_wl_borders[[0, -1]], 0)
        # s1 is the 0-based pixel that contains the lower limit
        # s2 is the 0-based pixel that contains the upper limit
        s1, s2 = ss_vals
        s1 = utils.coor_to_pix_1d(s1)
        s2 = utils.coor_to_pix_1d(s2)
        lower = max(0, min(s1, npix - 1))
        upper = max(0, min(s2, npix - 1))

        if lower > upper:
            warnings.warn('lower limit is > upper limit', RuntimeWarning)

        # We need a monotonic interpolator
        # linear would work, we use a cubic interpolator
        interpolator = SteffenInterpolator(
            old_wl_borders,
            accum_flux[idx],
            extrapolate='border'
        )

        if lower + span > upper - span:
            warnings.warn('lower limit + span is > upper limit - span', RuntimeWarning)

        fl_borders = interpolator(new_wl_borders)
        rss_resampled[idx] = fl_borders[1:] - fl_borders[:-1]
        # Expand the border to remove `span` pixels
        # in both sides, to avoid high variance
        rss_resampled[idx, lower:lower + span] = fill
        rss_resampled[idx, upper + 1 - span:upper + 1] = fill
        limits.append((fibid, (lower + span, upper - span)))

    return rss_resampled, limits


def pixel_borders(arr):
    import numina.array.wavecalib.resample as W
    return W.map_borders(arr)


if __name__ == '__main__':
    import numina.types.structured as stru
    import logging

    logging.basicConfig(level=logging.DEBUG)

    rssname = "/home/spr/Documentos/Congresos/ria_megara/M15_LCB_HR-R/obsid3_results/reduced_rss.fits"
    rss = fits.open(rssname)
    wlname = "/home/spr/Documentos/Congresos/ria_megara/M15_LCB_HR-R/obsid2_results/master_wlcalib.json"
    solutionwl = stru.open(wlname)

    rss2 = calibrate_wl_rss_megara(rss, solutionwl)
    rss2.writeto('result_span0b.fits', overwrite=True)

    rss2 = calibrate_wl_rss_megara(rss, solutionwl, span=2)
    rss2.writeto('result_span2b.fits', overwrite=True)