"""The recipes of the combined modes, not implemented yet"""

import pytest

from megaradrp.loader import load_drp
from megaradrp.recipes.combined import extinctionstar, sensstar


@pytest.mark.parametrize("module", [sensstar, extinctionstar])
def test_not_implemented(module):
    with pytest.raises(NotImplementedError, match="not implemented yet"):
        module.Recipe().run(None)


@pytest.mark.parametrize("mode", ["MegaraSensitivityStar", "MegaraExtinctionStar"])
def test_not_in_drp(mode):
    drp = load_drp()
    assert mode not in drp.modes
    assert mode not in drp.pipelines["default"].recipes
