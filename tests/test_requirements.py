import pytest
from numina.core import Requirement

from megaradrp.requirements import EXTRACTION_OFFSET_AUTO, DiffuseLightRequirement, ExtractionOffsetParameter


def test_requires_df_l():
    req = DiffuseLightRequirement()
    assert isinstance(req, Requirement)


@pytest.mark.parametrize(
    "value, expected",
    [("auto", "auto"), (" AUTO ", "auto"), ([1.0, 2e-5], [1.0, 2e-5]), (1.5, [1.5])],
)
def test_extraction_offset_parameter(value, expected):
    param = ExtractionOffsetParameter()
    converted = param.convert(value)
    assert converted == expected
    assert param.validate(converted)


def test_extraction_offset_parameter_default_and_invalid():
    param = ExtractionOffsetParameter()
    assert param.default_value() == [0.0]
    with pytest.raises(ValueError, match=EXTRACTION_OFFSET_AUTO):
        param.convert("yes")


def test_extraction_offset_load():
    """The values of the observation result are loaded with the type"""
    import numina.store

    param = ExtractionOffsetParameter()
    assert numina.store.load(param.type, "auto") == "auto"
    assert numina.store.load(param.type, [1.0, 2.0]) == [1.0, 2.0]
