import types

from astropy.io import fits
import numpy
import pytest
from numina.array.crmasks.valid_parameters import VALID_COMBINATIONS
from numina.core import DataFrame

import megaradrp.recipes.calibration.crdetect as crdetect
from megaradrp.recipes.calibration.crdetect import CRMasksRecipe

NIMAGES = 3


def make_frame(idx):
    hdu = fits.PrimaryHDU(numpy.full((10, 20), 100.0 + idx, dtype="float32"))
    hdu.header["UUID"] = f"00000000-0000-0000-0000-00000000000{idx}"
    return DataFrame(frame=fits.HDUList([hdu]))


def overscan_trimming(hdul):
    hdul[0].header.add_history("overscan and trimming")
    return hdul


def bias_dark_gain(hdul):
    hdul[0].data = hdul[0].data * 2
    hdul[0].header.add_history("bias, dark and gain")
    return hdul


@pytest.fixture
def calls(monkeypatch):
    """Replace the reduction flows and the computation of the masks, record the calls"""
    calls = {"compute": [], "apply": []}
    monkeypatch.setattr(CRMasksRecipe, "init_filters", lambda self, rinput, conf: [overscan_trimming, bias_dark_gain])

    def compute_crmasks(list_arrays, **kwds):
        calls["compute"].append((list_arrays, kwds))
        shape = list_arrays[0].shape
        hdus = [fits.PrimaryHDU()]
        # extension i has i + 1 masked pixels
        names = ["MEDIANCR", "MEANCRT"] + [f"CRMASK{i + 1}" for i in range(len(list_arrays))]
        for i, name in enumerate(names):
            mask = numpy.zeros(shape, dtype="uint8")
            mask.flat[: i + 1] = 1
            hdus.append(fits.ImageHDU(mask, name=name))
        return fits.HDUList(hdus)

    def apply_crmasks(list_arrays, hdul_masks, combination, **kwds):
        calls["apply"].append((combination, kwds))
        return numpy.sum(list_arrays, axis=0), None, None

    monkeypatch.setattr(crdetect, "compute_crmasks", compute_crmasks)
    monkeypatch.setattr(crdetect, "apply_crmasks", apply_crmasks)
    return calls


def make_input(**kwds):
    """The default values of the parameters of the recipe"""
    values = {key: req.default for key, req in CRMasksRecipe.requirements().items()}
    values["obresult"] = types.SimpleNamespace(frames=[make_frame(idx) for idx in range(NIMAGES)], configuration=None)
    values.update(kwds)
    return types.SimpleNamespace(**values)


def test_crdetect(calls):
    result = CRMasksRecipe().run(make_input(rnoise=4.0, crmethod="mm_lacosmic", interactive=False))
    ((arrays, kwds),) = calls["compute"]
    # the images are reduced with the two flows
    assert len(arrays) == NIMAGES
    for idx, array in enumerate(arrays):
        assert numpy.all(array == 2 * (100.0 + idx))
    # the images are in electrons, without bias
    assert (kwds["gain"], kwds["bias"], kwds["rnoise"]) == (1.0, 0.0, 4.0)
    assert kwds["crmethod"] == "mm_lacosmic"
    assert kwds["interactive"] is False
    # all the other parameters are passed with their names
    for key in ["la_sigclip", "pc_rlim", "dc_mask", "nn_model", "mm_seed"]:
        assert kwds[key] == CRMasksRecipe.requirements()[key].default
    assert calls["apply"] == []

    hdul = result.crmasks.open()
    history = list(hdul[0].header["HISTORY"])
    for idx in range(NIMAGES):
        assert f"Image 00000000-0000-0000-0000-00000000000{idx}" in history
    assert history.count("bias, dark and gain") == NIMAGES
    assert "Extension MEDIANCR: 1 masked pixels" in history
    assert f"Extension CRMASK{NIMAGES}: {NIMAGES + 2} masked pixels" in history
    assert hdul[0].header["NUMRNAM"] == "CRMasksRecipe"


def test_crdetect_intermediate(calls, tmp_path, monkeypatch):
    """The preprocessed and postprocessed images are saved"""
    monkeypatch.chdir(tmp_path)
    recipe = CRMasksRecipe(intermediate_results=True)
    recipe.run(make_input(save_preprocessed=True, save_postprocessed=True))
    for idx in range(NIMAGES):
        with fits.open(tmp_path / f"preprocessed_{idx + 1}.fits") as hdul:
            assert numpy.all(hdul[0].data == 2 * (100.0 + idx))
    with fits.open(tmp_path / "preprocessed_median.fits") as hdul:
        assert numpy.all(hdul[0].data == 202.0)
        assert hdul[0].data.dtype.name == "float32"
    assert [combination for combination, _ in calls["apply"]] == VALID_COMBINATIONS
    for combination in VALID_COMBINATIONS:
        assert (tmp_path / f"postprocessed_{combination}.fits").exists()


def test_crdetect_invalid_flow(calls, monkeypatch):
    monkeypatch.setattr(CRMasksRecipe, "init_filters", lambda self, rinput, conf: [overscan_trimming])
    with pytest.raises(ValueError, match="Invalid reduction flow"):
        CRMasksRecipe().run(make_input())
