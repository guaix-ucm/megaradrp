from types import SimpleNamespace

import astropy.io.fits as fits
import jsonschema
import pytest

from megaradrp.datatype import MegaraDataType
from megaradrp.validators import ExtChecker, check_as_datatype, definition_validator


def create_hdulist(*cards_list):
    # The checkers only read the headers of the HDUs
    return [SimpleNamespace(header=fits.Header(list(cards.items()))) for cards in cards_list]


def bias_cards(**kwargs):
    cards = {
        "NAXIS": 2,
        "NAXIS1": 4196,
        "NAXIS2": 4212,
        "INSTRUME": "MEGARA",
        "OBSMODE": "MegaraBiasImage",
        "EXPTIME": 0.0,
        "DARKTIME": 0.0,
    }
    cards.update(kwargs)
    return cards


def test_definitions_exist():
    definitions = check_as_datatype.validator_image.schema["definitions"]
    for checker in check_as_datatype._megara_checkers.values():
        if isinstance(checker, ExtChecker):
            for name in checker.definitions:
                assert name in definitions


def test_check_bias():
    checker = check_as_datatype(MegaraDataType.IMAGE_BIAS)
    assert checker(create_hdulist(bias_cards()))


@pytest.mark.parametrize(
    "cards",
    [
        {"EXPTIME": 10.0},
        {"OBSMODE": "MegaraDarkImage"},
        {"IMAGETYP": "IMAGE_DARK"},
        {"NAXIS1": 4096},
    ],
)
def test_check_bias_invalid(cards):
    checker = check_as_datatype(MegaraDataType.IMAGE_BIAS)
    with pytest.raises(jsonschema.exceptions.ValidationError):
        checker(create_hdulist(bias_cards(**cards)))


def test_check_bias_number_of_hdus():
    checker = check_as_datatype(MegaraDataType.IMAGE_BIAS)
    with pytest.raises(ValueError, match="number of HDUs"):
        checker(create_hdulist(bias_cards(), {"EXTNAME": "FIBERS"}))


def test_definition_validator_resolves_references():
    # rss_hdu_values refers to rss_naxis2
    validator = definition_validator(check_as_datatype.validator_image, "rss_hdu_values")
    validator.validate({"NAXIS": 2, "NAXIS1": 4300, "NAXIS2": 623, "INSMODE": "LCB"})
    with pytest.raises(jsonschema.exceptions.ValidationError):
        validator.validate({"NAXIS": 2, "NAXIS1": 4300, "NAXIS2": 644, "INSMODE": "LCB"})
