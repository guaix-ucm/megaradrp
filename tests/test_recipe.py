import types

import numpy
import pytest
from astropy.io import fits
from numina.core import BaseRecipe, DataFrame
from numina.core.recipeinout import TagsChecker
from numina.exceptions import ValidationError
from numina.types.qc import QC

import megaradrp.requirements as reqs
from megaradrp.core.recipe import MegaraBaseRecipe, MegaraRecipeInput
from megaradrp.testing.common import make_obresult, raw_header
from megaradrp.testing.create_tracemap import create_test_tracemap2
from megaradrp.testing.create_traceoffset import create_offset_image, create_offset_tracemap


def test_base_recipe():
    version = "1.0.1"
    obj = MegaraBaseRecipe(version=version)
    assert isinstance(obj, BaseRecipe)
    assert obj.__version__ == version


class TagsRecipe(MegaraBaseRecipe):
    """A recipe with requirements with tags, one of them a MultiType and one optional"""

    master_traces = reqs.MasterTraceMapRequirement()
    master_apertures = reqs.MasterAperturesRequirement()
    master_slitflat = reqs.MasterSlitFlatRequirement()

    def run(self, rinput):
        return self.create_result()


def make_tracemap(vph):
    tracemap = create_test_tracemap2()
    tracemap.uuid = "410f6ea0-c3df-481c-9820-5cf9e0ed3d91"
    tracemap.tags = {"insmode": "LCB", "vph": vph}
    return tracemap


@pytest.fixture(scope="module")
def obresult():
    obresult = make_obresult("MegaraTraceMap", header=raw_header(vph="LR-R", insmode="LCB"))
    obresult.tags = {"vph": "LR-R", "insmode": "LCB"}
    return obresult


def test_validate_input(obresult):
    """The tags of the requirements are those of the observation result"""
    recipe = TagsRecipe()
    rinput = recipe.create_input(
        obresult=obresult, master_traces=make_tracemap("LR-R"), master_apertures=make_tracemap("LR-R")
    )
    recipe.validate_input(rinput)


@pytest.mark.parametrize(
    "traces_vph, apertures_vph, invalid",
    [
        ("LR-B", "LR-R", "master_traces"),
        ("LR-R", "LR-B", "master_apertures"),
        # the first one is reported
        ("LR-B", "LR-B", "master_traces"),
    ],
)
def test_validate_input_invalid_tags(obresult, traces_vph, apertures_vph, invalid):
    recipe = TagsRecipe()
    rinput = recipe.create_input(
        obresult=obresult, master_traces=make_tracemap(traces_vph), master_apertures=make_tracemap(apertures_vph)
    )
    msg = f"requirement {invalid}: the value of {invalid} has vph=LR-B, the observing block has vph=LR-R"
    with pytest.raises(ValidationError, match=msg):
        recipe.validate_input(rinput)


def test_recipes_check_tags():
    """The inputs of the recipes check the tags of their products"""
    assert issubclass(TagsRecipe.RecipeInput, MegaraRecipeInput)
    assert any(isinstance(checker, TagsChecker) for checker in TagsRecipe.RecipeInput.__checkers__)


def test_validate_input_types(obresult):
    """The values are validated by their types before the tags"""
    recipe = TagsRecipe()
    tracemap = make_tracemap("LR-R")
    tracemap.uuid = "123"
    rinput = recipe.create_input(obresult=obresult, master_traces=tracemap, master_apertures=make_tracemap("LR-R"))
    with pytest.raises(ValidationError, match="master_traces: uuid"):
        recipe.validate_input(rinput)


def test_run_qc():
    recipe = TagsRecipe()
    result = recipe.run_qc(None, types.SimpleNamespace(qc=QC.UNKNOWN))
    assert result.qc == QC.GOOD


def test_recipe_logger():
    """The recipes log with the logger of their module"""
    from megaradrp.recipes.calibration.arc import ArcCalibrationRecipe

    assert ArcCalibrationRecipe.logger.name == "megaradrp.recipes.calibration.arc"
    assert TagsRecipe.logger.name == TagsRecipe.__module__


class OffsetRecipe(MegaraBaseRecipe):
    extraction_offset = reqs.ExtractionOffsetParameter()

    def run(self, rinput):
        return self.create_result()


def test_extraction_offset_input(obresult):
    recipe = OffsetRecipe()
    assert recipe.create_input(obresult=obresult).extraction_offset == [0.0]
    assert recipe.create_input(obresult=obresult, extraction_offset="auto").extraction_offset == "auto"
    assert recipe.create_input(obresult=obresult, extraction_offset=[1.0, 2e-4]).extraction_offset == [1.0, 2e-4]


def offset_frames(tracemap, offsets, **kwds):
    frames = []
    for i, c0 in enumerate(offsets):
        header = fits.Header({"RDNOISE1": 3.4, "UUID": f"00000000-0000-4000-8000-00000000000{i}"})
        image = create_offset_image(tracemap, c0, seed=i, **kwds)
        frames.append(DataFrame(frame=fits.HDUList([fits.PrimaryHDU(image, header)])))
    return frames


def test_resolve_extraction_offset():
    """With "auto", the mean of the offsets measured in each frame"""
    tracemap = create_offset_tracemap()
    frames = offset_frames(tracemap, [1.0, 1.1])
    rinput = types.SimpleNamespace(extraction_offset="auto", obresult=types.SimpleNamespace(frames=frames))
    coef, history = MegaraBaseRecipe().resolve_extraction_offset(rinput, tracemap, flow=[])
    y_mid = 165.0
    assert coef[0] + coef[1] * y_mid == pytest.approx(1.05, abs=0.03)
    assert len(history) == 3
    assert history[0].startswith("Extraction offset of image 1: [1.0")
    assert history[-1].startswith("Extraction offset measured:")


def test_resolve_extraction_offset_given():
    rinput = types.SimpleNamespace(extraction_offset=[0.5, 1e-4])
    assert MegaraBaseRecipe().resolve_extraction_offset(rinput, None) == ([0.5, 1e-4], [])


def test_resolve_extraction_offset_no_signal(caplog):
    """Without a reliable measurement, the offset is 0"""
    tracemap = create_offset_tracemap()
    frames = offset_frames(tracemap, [1.0], flux=numpy.zeros(40))
    rinput = types.SimpleNamespace(extraction_offset="auto", obresult=types.SimpleNamespace(frames=frames))
    coef, history = MegaraBaseRecipe().resolve_extraction_offset(rinput, tracemap, flow=[])
    assert coef == [0.0]
    assert "using 0" in history[-1]
    assert "no reliable measurement" in caplog.text
