import astropy.io.fits as fits
import numpy
import numina.core.config as cfg
import numina.core.pipeline as ppl
import pytest

import megaradrp.datamodel as DM
import megaradrp.loader as loader
from megaradrp.datatype import DataOrigin
from megaradrp.loader import load_drp
from megaradrp.testing.common import raw_header


def test_loader():

    drpm = load_drp()
    assert isinstance(drpm, ppl.InstrumentDRP)
    assert drpm.name == "MEGARA"


@pytest.fixture
def megara_fits(tmp_path):
    hdu = fits.PrimaryHDU(numpy.zeros((10, 10), dtype="uint16"), header=raw_header(vph="LR-R", insmode="LCB"))
    hdu.header["UUID"] = "410f6ea0-c3df-481c-9820-5cf9e0ed3d91"
    hdu.header["IMAGETYP"] = "BIAS"
    path = tmp_path / "megara.fits"
    hdu.writeto(path)
    return path


def test_is_fits_megara(megara_fits, tmp_path):
    assert loader.is_fits_megara(str(megara_fits))
    other = tmp_path / "other.fits"
    fits.writeto(other, numpy.zeros((2, 2)))
    assert not loader.is_fits_megara(str(other))
    assert not loader.is_fits_megara(str(tmp_path / "other.json"))


def test_describe_fits_megara(megara_fits):
    """The images of MEGARA are described by the function of the DRP"""
    description = cfg.describe(str(megara_fits))
    assert description["instrument"] == "MEGARA"
    assert description["uuid"] == "410f6ea0-c3df-481c-9820-5cf9e0ed3d91"
    assert description["insconf"] == raw_header()["INSCONF"]
    assert description["origin"] == DataOrigin.OBSERVED


def test_check_obj_megara(monkeypatch):
    """The checker of the DRP is the checker of the datamodel"""
    calls = []
    monkeypatch.setattr(DM, "check_obj_megara", lambda obj, astype, level: calls.append((obj, astype, level)))
    drpm = load_drp()
    drpm.checker("obj", astype="type", level=2)
    assert calls == [("obj", "type", 2)]
