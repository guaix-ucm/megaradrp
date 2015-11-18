import inspect
import sys
from megaradrp.recipes.calibration.cBase import MegaraBaseRecipe
from numina.core import BaseRecipeAutoQC

def test_cBase(benchmark):
    version = "1.0.1"
    obj = MegaraBaseRecipe(version)
    assert isinstance(obj,BaseRecipeAutoQC)
    assert obj.__version__ == version

