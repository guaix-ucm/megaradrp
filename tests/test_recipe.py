import types

import pytest
from numina.core import BaseRecipe
from numina.core.recipeinout import TagsChecker
from numina.exceptions import ValidationError
from numina.types.qc import QC

import megaradrp.requirements as reqs
from megaradrp.core.recipe import MegaraBaseRecipe, MegaraRecipeInput
from megaradrp.testing.common import make_obresult, raw_header
from megaradrp.testing.create_tracemap import create_test_tracemap2


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
