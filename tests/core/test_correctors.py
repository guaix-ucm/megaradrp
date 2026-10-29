import types

import astropy.io.fits as fits
import numpy
import numina.processing as proc
import pytest
from numina.core import DataFrame
from numina.util.node import IdNode

import megaradrp.core.correctors as cor
from megaradrp.datamodel import MegaraDataModel
from megaradrp.processing.diffuselight import DiffuseLightCorrector
from megaradrp.processing.slitflat import SlitFlatCorrector
from megaradrp.processing.trimover import GainCorrector, OverscanCorrector, TrimImage
from megaradrp.testing.common import make_obresult

UUID = "410f6ea0-c3df-481c-9820-5cf9e0ed3d91"


def calibration(value=2.0):
    hdu = fits.PrimaryHDU(numpy.full((10, 20), value, dtype="float32"))
    hdu.header["UUID"] = UUID
    return DataFrame(frame=fits.HDUList([hdu]))


# key of the requirement, function, class of the corrector
CORRECTORS = [
    ("master_bpm", cor.get_corrector_bpm, proc.BadPixelCorrector),
    ("master_bias", cor.get_corrector_bias, proc.BiasCorrector),
    ("master_dark", cor.get_corrector_dark, proc.DarkCorrector),
    ("master_slitflat", cor.get_corrector_slit_flat, SlitFlatCorrector),
    ("diffuse_light_image", cor.get_corrector_diffuse_light, DiffuseLightCorrector),
]


@pytest.mark.parametrize("key, getter, corrector_class", CORRECTORS)
def test_get_corrector(key, getter, corrector_class):
    """The corrector of the calibration, with its id"""
    rinput = types.SimpleNamespace(**{key: calibration()})
    meta = {key: {"uuid": UUID}}
    corrector = getter(rinput, meta, None, MegaraDataModel())
    assert isinstance(corrector, corrector_class)
    assert corrector.calibid == UUID


@pytest.mark.parametrize("key, getter, corrector_class", CORRECTORS)
def test_get_corrector_missing(key, getter, corrector_class):
    """Without the calibration, the corrector does nothing"""
    rinput = types.SimpleNamespace(**{key: None})
    corrector = getter(rinput, {}, None, MegaraDataModel())
    assert isinstance(corrector, IdNode)


def test_get_corrector_bpm_window():
    """The bad pixels are interpolated along the rows"""
    rinput = types.SimpleNamespace(master_bpm=calibration(0.0))
    corrector = cor.get_corrector_bpm(rinput, {}, None, MegaraDataModel())
    assert (corrector.hwin, corrector.wwin) == (0, 3)


def test_get_corrector_super_deprecated():
    rinput = types.SimpleNamespace(master_bias=calibration())
    with pytest.warns(DeprecationWarning, match="use numina.processing.getters.load_corrector"):
        corrector = cor.get_corrector_super(rinput, {}, "master_bias", proc.BiasCorrector, MegaraDataModel())
    assert isinstance(corrector, proc.BiasCorrector)
    assert corrector.calibid == UUID


def test_get_corrector_slit_flat_applied():
    """The image is divided by the slit flat"""
    rinput = types.SimpleNamespace(master_slitflat=calibration(2.0))
    corrector = cor.get_corrector_slit_flat(rinput, {"master_slitflat": {}}, None, MegaraDataModel())
    img = fits.HDUList([fits.PrimaryHDU(numpy.full((10, 20), 10.0, dtype="float32"))])
    img[0].header["UUID"] = "11111111-2222-3333-4444-555555555555"
    result = corrector(img)
    assert numpy.allclose(result[0].data, 5.0)


@pytest.fixture(scope="module")
def instrument():
    return make_obresult("MegaraBiasImage").configuration


@pytest.mark.parametrize(
    "getter, corrector_class",
    [
        (cor.get_corrector_overscan, OverscanCorrector),
        (cor.get_corrector_trimming, TrimImage),
        (cor.get_corrector_gain, GainCorrector),
    ],
)
def test_get_corrector_detector(instrument, getter, corrector_class):
    """The correctors of the detector have the id of the detector"""
    corrector = getter(None, {}, instrument, MegaraDataModel())
    assert isinstance(corrector, corrector_class)
    assert corrector.calibid == str(instrument.get_device("detector").origin.uuid)
