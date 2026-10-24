import pytest

from astropy.io import fits
import numpy

from megaradrp.datamodel import MegaraDataModel
from megaradrp.processing.fluxcalib import FluxCalibration, update_flux_limits, WAVLIM_KEYS, PIXLIM_KEYS
from megaradrp.testing.create_wcs import generate_wcs


@pytest.mark.parametrize("wcsl", [None, generate_wcs()])
@pytest.mark.parametrize("ref", [0, 1])
def test_update_flux_limits2(wcsl, ref):

    header = {}

    values = numpy.array([1, 1000, 2, 900, 100, 800])

    off = (ref + 1) % 2

    if wcsl is None:
        results = values + off
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


def make_sensitivity(nwl=50):
    sens = numpy.zeros(nwl, dtype="float32")
    sens[5:45] = 4.0
    hdu = fits.PrimaryHDU(sens)
    hdr = hdu.header
    hdr["UUID"] = "11111111-2222-3333-4444-555555555555"
    hdr["TUNIT"] = "Jy"
    # 1-based limits
    for key, val in zip(PIXLIM_KEYS, [6, 45, 6, 45, 10, 40]):
        hdr[key] = val
    return fits.HDUList([hdu])


def make_rss(nwl=50):
    hdu = fits.PrimaryHDU(numpy.full((3, nwl), 40.0, dtype="float32"))
    hdr = hdu.header
    hdr["UUID"] = "410f6ea0-c3df-481c-9820-5cf9e0ed3d91"
    hdr["EXPTIME"] = 2.0
    hdr["CTYPE1"] = "AWAV"
    hdr["CRPIX1"] = 1
    hdr["CRVAL1"] = 4000.0
    hdr["CDELT1"] = 2.0
    hdr["CUNIT1"] = "Angstrom"
    return fits.HDUList([hdu])


def test_flux_calibration():
    sens = make_sensitivity()
    corrector = FluxCalibration(sens, datamodel=MegaraDataModel())
    img = corrector(make_rss())
    data = img[0].data
    hdr = img[0].header

    # counts / exptime / sensitivity, in the region with at least one fiber
    assert numpy.allclose(data[:, 5:45], 40.0 / 2.0 / 4.0)
    assert numpy.all(data[:, :5] == 0)
    assert numpy.all(data[:, 45:] == 0)
    assert hdr["BUNIT"] == "Jy"
    assert hdr["NUM-FLUX"] == "11111111-2222-3333-4444-555555555555"
    # the limits of the sensitivity, 1-based, and their wavelengths
    for key, val in zip(PIXLIM_KEYS, [6, 45, 6, 45, 10, 40]):
        assert hdr[key] == val
        assert hdr[key.replace("PIX", "WAV")] == pytest.approx(4000.0 + (val - 1) * 2.0)
