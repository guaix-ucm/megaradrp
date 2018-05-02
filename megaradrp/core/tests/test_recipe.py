

from numina.core import BaseRecipe

from megaradrp.core.recipe import MegaraBaseRecipe


def test_base_recipe():
    version = "1.0.1"
    obj = MegaraBaseRecipe(version=version)
    assert isinstance(obj, BaseRecipe)
    assert obj.__version__ == version