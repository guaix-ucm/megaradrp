import pytest

import numpy
import astropy.wcs

from ..fluxcalib import update_flux_limits, WAVLIM_KEYS, PIXLIM_KEYS


def generate_wcs():
    wcsl = astropy.wcs.WCS(naxis=2)

    crpix = 3.0
    crval = 5300.0
    cdelt = 0.1

    wcsl.wcs.crpix = [crpix, 0.0]
    wcsl.wcs.crval = [crval, 0.0]
    wcsl.wcs.cdelt = [cdelt, 1.0]
    return wcsl


@pytest.mark.parametrize("wcsl", [None, generate_wcs()])
@pytest.mark.parametrize("ref", [0, 1])
def test_update_flux_limits2(wcsl, ref):

    header = {}

    values = numpy.array([1, 1000, 2, 900, 100, 800])

    off = (ref + 1) % 2

    if wcsl is None:
        results = (values + off)
    else:
        lm = numpy.array([values, numpy.zeros_like(values)])
        wavelen_ = wcsl.all_pix2world(lm.T, ref)
        results = wavelen_[:, 0]

    fluxlimits = dict(zip(PIXLIM_KEYS, values))
    hdr = update_flux_limits(header, fluxlimits, wcs=wcsl, ref=ref)

    for key in WAVLIM_KEYS:
        assert key in hdr

    for key in PIXLIM_KEYS:
        assert key in hdr

    for key, val in zip(PIXLIM_KEYS, values):
        assert hdr[key][0] == val + off

    for key, res in zip(WAVLIM_KEYS, results):
        assert hdr[key][0] == res


def test_update_flux_invalid():

    header = {}

    values = [1, 1000, 2, 900, 100, 800]
    fluxlimits = dict(zip(PIXLIM_KEYS, values))

    ref = 2
    with pytest.raises(ValueError):
        update_flux_limits(header, fluxlimits, ref=ref)
