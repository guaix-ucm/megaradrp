import pytest
import astropy.io.fits as fits
import numpy
import numpy.polynomial.polynomial as nppol
from numina.array.wavecalib.solutionarc import CrLinear, SolutionArcCalibration

from megaradrp.instrument import WLCALIB_PARAMS
from megaradrp.products.wavecalibration import FiberSolutionArcCalibration, WavelengthCalibration
from megaradrp.testing.create_header import create_spec_header
from megaradrp.testing.create_image import create_rss
from megaradrp.processing.wavecalibration import (
    SimpleWcs1D,
    WavelengthCalibrator,
    calibrate_wl_rss,
    calibrate_wl_rss_megara,
    header_add_barycentric_correction,
    resample_rss_flux,
)


def test_add_barycentric_missing1():
    hdr = create_spec_header()
    del hdr["RADEG"]

    with pytest.raises(KeyError):
        header_add_barycentric_correction(hdr, key="b")


def test_add_barycentric_missing2():
    hdr = create_spec_header()
    del hdr["DATE-OBS"]

    with pytest.raises(KeyError):
        header_add_barycentric_correction(hdr, key="b")


def test_add_barycentric_missing3():
    hdr = fits.Header()
    hdr["DATE-OBS"] = "2017-08-23T21:38:30.55"
    # GTC
    # hdr["OBSGEO-X"] = 5327285.0921
    # hdr["OBSGEO-Y"] = -1718777.1125
    # hdr["OBSGEO-Z"] = 3051786.7327
    hdr["OBSGEO-B"] = +28.76060
    hdr["OBSGEO-L"] = -17.88160
    hdr["OBSGEO-H"] = 2326
    hdr["RADEG"] = 285.481037748898
    hdr["DECDEG"] = 42.4882140636786

    with pytest.raises(TypeError):
        header_add_barycentric_correction(hdr, key="b")


def make_solution(fibids, coeff, offset=0.0):
    """A wavelength calibration with the same polynomial in all the fibers"""
    wlcalib = WavelengthCalibration(instrument="MEGARA")
    for fibid in fibids:
        sol = SolutionArcCalibration([], list(coeff), 0.0, CrLinear(0, 0, 0, 0, 0))
        wlcalib.contents.append(FiberSolutionArcCalibration(fibid, sol))
    wlcalib.global_offset = nppol.Polynomial([offset])
    wlcalib.error_fitting = []
    wlcalib.missing_fibers = []
    return wlcalib


# 1-based pixel x of the input has wavelength 4000 + x
COEFF = [4000.0, 1.0]


@pytest.fixture
def raw_rss():
    rng = numpy.random.default_rng(0)
    return rng.uniform(1, 2, size=(3, 200))


def test_resample_rss_flux(raw_rss):
    """The flux is conserved when the target range contains the input"""
    solution = make_solution([1, 2, 3], COEFF)
    target = SimpleWcs1D(crval=3990.0, crpix=1.0, cdelt=0.5).create_internal_wcs_()
    result, limits = resample_rss_flux(raw_rss, solution, 500, target)

    assert result.shape == (3, 500)
    assert numpy.allclose(result.sum(axis=1), raw_rss.sum(axis=1))
    # the input covers from 4000.5 to 4200.5, the pixels containing the limits
    for fibid, (lower, upper) in limits:
        assert (lower, upper) == (21, 421)
        assert numpy.all(result[fibid - 1, :lower] == 0)
        assert numpy.all(result[fibid - 1, upper + 1 :] == 0)


def test_resample_rss_flux_span(raw_rss):
    """The span pixels at both sides are filled and removed from the limits"""
    solution = make_solution([1, 2, 3], COEFF)
    target = SimpleWcs1D(crval=3990.0, crpix=1.0, cdelt=0.5).create_internal_wcs_()
    result, limits = resample_rss_flux(raw_rss, solution, 500, target, span=2, fill=-1)

    assert limits[0] == (1, (23, 419))
    assert numpy.all(result[:, 21:23] == -1)
    assert numpy.all(result[:, 420:422] == -1)


def test_resample_rss_flux_offset(raw_rss):
    """The global offset is subtracted from the wavelength of the solution"""
    target = SimpleWcs1D(crval=3990.0, crpix=1.0, cdelt=0.5).create_internal_wcs_()
    _, limits = resample_rss_flux(raw_rss, make_solution([1], COEFF, offset=3.0), 500, target)
    # 3 Angstrom, 6 pixels to the blue
    assert limits[0] == (1, (15, 415))


@pytest.mark.xfail(reason="the global offset is subtracted in the solution, again in each call", strict=True)
def test_resample_rss_flux_twice(raw_rss):
    """Applying a solution twice gives the same result"""
    solution = make_solution([1, 2, 3], COEFF, offset=3.0)
    target = SimpleWcs1D(crval=3990.0, crpix=1.0, cdelt=0.5).create_internal_wcs_()
    result1, limits1 = resample_rss_flux(raw_rss, solution, 500, target)
    result2, limits2 = resample_rss_flux(raw_rss, solution, 500, target)
    assert limits1 == limits2
    assert numpy.allclose(result1, result2)


def make_raw_rss(raw_rss):
    data = numpy.zeros((623, 200), dtype="float32")
    data[:3] = raw_rss
    rss = create_rss(data, numpy.zeros_like(data))
    # without WLMAP, it is created by the calibration
    del rss["WLMAP"]
    rss[0].header["UUID"] = "410f6ea0-c3df-481c-9820-5cf9e0ed3d91"
    rss[0].header["INSMODE"] = "LCB"
    rss[0].header["VPH"] = "LR-B"
    return rss


def test_calibrate_wl_rss(raw_rss):
    rss = make_raw_rss(raw_rss)
    solution = make_solution([1, 2, 3], COEFF)
    solution.uuid = "11111111-2222-3333-4444-555555555555"
    solution.error_fitting = [2]
    solution.missing_fibers = [3]
    target = SimpleWcs1D(crval=3990.0, crpix=1.0, cdelt=0.5)
    result = calibrate_wl_rss(rss, solution, 500, target)

    # a copy, the input is not modified
    assert result is not rss
    assert rss[0].data.shape == (623, 200)
    hdr = result[0].header
    assert result[0].data.shape == (623, 500)
    assert result[0].data.dtype == numpy.float32
    assert (hdr["CTYPE1"], hdr["CRVAL1"], hdr["CDELT1"], hdr["CRPIX1"]) == ("AWAV", 3990.0, 0.5, 1.0)
    assert hdr["CUNIT1"] == "Angstrom"
    assert hdr["NUM-WAV"] == "uuid:11111111-2222-3333-4444-555555555555"
    assert hdr["UUID"] != "410f6ea0-c3df-481c-9820-5cf9e0ed3d91"
    # limits of the spectral coverage, 1-based, and the map of valid pixels
    fibers = result["FIBERS"].header
    assert (fibers["FIB001W1"], fibers["FIB001W2"]) == (22, 422)
    wlmap = result["WLMAP"].data
    assert numpy.all(wlmap[0, 21:422] == 1)
    assert wlmap[0].sum() == 401
    assert wlmap[3:].sum() == 0
    # the fibers with errors in the calibration are invalid
    assert fibers["FIB002_V"] is False
    assert fibers["FIB003_V"] is False
    assert "FIB001_V" not in fibers


def test_calibrate_wl_rss_megara(raw_rss):
    """The target WCS is the one of the mode and VPH"""
    rss = make_raw_rss(raw_rss)
    solution = make_solution([1], [4280.0, 1.0])
    result = calibrate_wl_rss_megara(rss, solution)
    params = WLCALIB_PARAMS["LCB"]["LR-B"]
    assert result[0].data.shape == (623, params["npix"])
    assert result[0].header["CRVAL1"] == params["crval"]
    assert result[0].header["CDELT1"] == params["cdelt"]


def test_calibrate_wl_rss_megara_undefined(raw_rss):
    rss = make_raw_rss(raw_rss)
    rss[0].header["VPH"] = "XX-X"
    with pytest.raises(ValueError):
        calibrate_wl_rss_megara(rss, make_solution([1], COEFF))


def test_wavelength_calibrator(raw_rss):
    """The corrector works in place, removing 2 pixels at both sides"""
    rss = make_raw_rss(raw_rss)
    solution = make_solution([1], [4280.0, 1.0])
    result = WavelengthCalibrator(solution)(rss)
    assert result is rss
    assert "Resample span=2" in result[0].header["HISTORY"]
