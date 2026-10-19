"""The requirements and results of the twilight recipe, declared in the recipe"""

from megaradrp.recipes.calibration import twilight
from megaradrp.recipes.calibration.twilight import TwilightFiberFlatRecipe


def test_fields():
    rinput = TwilightFiberFlatRecipe.RecipeInput
    assert list(rinput.stored()) == [
        "obresult",
        "master_bias",
        "master_dark",
        "master_bpm",
        "master_slitflat",
        "master_apertures",
        "extraction_offset",
        "normalize_region",
        "continuum_region",
        "master_wlcalib",
        "master_fiberflat",
    ]
    assert rinput.stored()["master_apertures"].alias == "master_traces"
    assert list(TwilightFiberFlatRecipe.RecipeResult.stored()) == [
        "reduced_image",
        "reduced_rss",
        "master_twilightflat",
    ]


def test_old_names():
    """The classes are still in the module with their old names"""
    assert twilight.RecipeInput is TwilightFiberFlatRecipe.RecipeInput
    assert twilight.RecipeResult is TwilightFiberFlatRecipe.RecipeResult
