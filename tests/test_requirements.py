from numina.core import Requirement
from megaradrp.requirements import DiffuseLightRequirement


def test_requires_df_l():
    req = DiffuseLightRequirement()
    assert isinstance(req, Requirement)
