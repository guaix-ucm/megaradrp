import pytest

from astropy.io import fits
import numpy

from megaradrp.datamodel import MegaraDataModel
from numina.exceptions import ValidationError

from megaradrp.processing.fluxcalib import FluxCalibration, update_flux_limits, WAVLIM_KEYS, PIXLIM_KEYS
from megaradrp.processing.fluxcalib import sensitivity_skip_windows, smooth_sensitivity, windows_check
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


NPIX = 2000
WL = 6000.0 + 0.5 * numpy.arange(NPIX)


def true_sensitivity():
    x = (WL - WL.mean()) / numpy.ptp(WL)
    return 1000.0 * (1.0 + 0.4 * x - 0.8 * x**2)


def raw_sensitivity(edge=100):
    """Raw sensitivity, 0 outside the region of valid flux, and its header"""
    raw = true_sensitivity()
    raw[:edge] = 0.0
    raw[-edge:] = 0.0
    header = fits.Header({"PIXLIMF1": edge + 1, "PIXLIMF2": NPIX - edge})
    return raw, header


def test_smooth_sensitivity_valid_region():
    """The zeros outside the region of valid flux do not bias the fit (#357); outside, straight lines"""
    raw, header = raw_sensitivity()
    smoothed, knots = smooth_sensitivity(WL, raw, header)
    region = slice(100, NPIX - 100)
    assert numpy.allclose(smoothed[region], true_sensitivity()[region], rtol=1e-6)
    # 5 fixed knots, equidistant in the region
    assert numpy.allclose(knots, numpy.linspace(WL[100], WL[NPIX - 101], 7)[1:-1])
    # outside the region, the value and the slope of the border
    assert numpy.allclose(numpy.diff(smoothed[:100], 2), 0.0, atol=1e-9)
    assert smoothed[99] == pytest.approx(smoothed[100], rel=1e-3)


def test_smooth_sensitivity_windows():
    """An absorption inside a skipped window does not change the fit"""
    raw, header = raw_sensitivity()
    band = (WL > 6400) & (WL < 6450)
    raw[band] *= 0.5
    smoothed, _ = smooth_sensitivity(WL, raw, header, skip_windows=[(6390.0, 6460.0)])
    region = slice(100, NPIX - 100)
    assert numpy.allclose(smoothed[region], true_sensitivity()[region], rtol=1e-6)
    biased, _ = smooth_sensitivity(WL, raw, header)
    assert numpy.max(numpy.abs(biased[region] / true_sensitivity()[region] - 1)) > 0.01


def test_smooth_sensitivity_reproducible():
    """With fixed knots, tiny changes of the wavelengths do not change the result (#389)"""
    raw, header = raw_sensitivity()
    raw = raw * (1 + 0.02 * numpy.sin(WL / 7.0))
    base, _ = smooth_sensitivity(WL, raw, header)
    for delta in (1e-12, -2e-12):
        other, _ = smooth_sensitivity(WL + delta, raw, header)
        assert numpy.allclose(other, base, rtol=1e-9)


def test_smooth_sensitivity_knots_list(caplog):
    """Knots outside the region of valid flux are not used"""
    raw, header = raw_sensitivity()
    smoothed, knots = smooth_sensitivity(WL, raw, header, knots=[6010.0, 6300.0, 6500.0, 6700.0])
    assert list(knots) == [6300.0, 6500.0, 6700.0]
    assert "knots outside the region of valid flux" in caplog.text


def test_smooth_sensitivity_too_many_windows(caplog):
    """Windows that leave too few pixels are ignored"""
    raw, header = raw_sensitivity()
    smoothed, _ = smooth_sensitivity(WL, raw, header, skip_windows=[(6000.0, 6900.0)])
    assert "they are ignored" in caplog.text
    region = slice(100, NPIX - 100)
    assert numpy.allclose(smoothed[region], true_sensitivity()[region], rtol=1e-6)


def test_smooth_sensitivity_undefined_with_windows(caplog):
    """If the windows leave knots without data, the fit is repeated without them"""
    raw, header = raw_sensitivity()
    knots = [6300.0, 6340.0, 6380.0, 6600.0]
    smoothed, _ = smooth_sensitivity(WL, raw, header, knots=knots, skip_windows=[(6290.0, 6390.0)])
    assert numpy.all(numpy.isfinite(smoothed))
    assert "repeating without them" in caplog.text


def test_windows_check():
    assert windows_check([[6000, 6100], (7000.0, 7100.5)]) == [[6000, 6100], (7000.0, 7100.5)]
    for bad in ([[6100, 6000]], [[6000]], [["a", 1]], [6000]):
        with pytest.raises(ValidationError):
            windows_check(bad)


def test_sensitivity_skip_windows():
    windows = sensitivity_skip_windows()
    assert all(w1 < w2 for w1, w2 in windows)
    # the O2 A band and H beta
    assert any(w1 < 7620 < w2 for w1, w2 in windows)
    assert any(w1 < 4861 < w2 for w1, w2 in windows)
    # the broad H2O band at the red end stays in the sensitivity
    assert not any(w1 < 9400 < w2 for w1, w2 in windows)
