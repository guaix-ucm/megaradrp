from types import SimpleNamespace

import astropy.io.fits as fits
import jsonschema
import pytest
from numina.exceptions import ValidationError

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


def lamp_hdulist(lamps_i=(), lamps_s=(), **kwargs):
    cards = {"NAXIS": 2, "NAXIS1": 4196, "NAXIS2": 4212}
    cards.update({f"LAMPI{idx}S": int(idx in lamps_i) for idx in range(1, 3)})
    cards.update({f"LAMPS{idx}S": int(idx in lamps_s) for idx in range(1, 6)})
    cards.update(kwargs)
    return create_hdulist(cards, {"EXTNAME": "FIBERS"})


@pytest.mark.parametrize(
    "datatype, lamps_i, lamps_s",
    [
        (MegaraDataType.IMAGE_FLAT, [1], []),
        (MegaraDataType.IMAGE_FLAT, [2], []),
        (MegaraDataType.IMAGE_COMP, [], [3, 5]),
        (MegaraDataType.IMAGE_TARGET, [], []),
    ],
)
def test_check_lamps(datatype, lamps_i, lamps_s):
    checker = check_as_datatype(datatype)
    assert checker(lamp_hdulist(lamps_i, lamps_s))


@pytest.mark.parametrize(
    "datatype, lamps_i, lamps_s, msg",
    [
        (MegaraDataType.IMAGE_FLAT, [], [], "all incandescent lamps are OFF"),
        (MegaraDataType.IMAGE_FLAT, [1], [2], "some comparison lamps are ON"),
        (MegaraDataType.IMAGE_COMP, [2], [1], "some incandescent lamps are ON"),
        (MegaraDataType.IMAGE_COMP, [], [], "all comparison lamps are OFF"),
        (MegaraDataType.IMAGE_TARGET, [1], [], "some incandescent lamps are ON"),
        (MegaraDataType.IMAGE_TARGET, [], [4], "some comparison lamps are ON"),
    ],
)
def test_check_lamps_invalid(datatype, lamps_i, lamps_s, msg):
    checker = check_as_datatype(datatype)
    with pytest.raises(ValidationError, match=msg):
        checker(lamp_hdulist(lamps_i, lamps_s))


def test_check_lamps_missing_keyword():
    hdulist = lamp_hdulist([1], [])
    del hdulist[0].header["LAMPS3S"]
    checker = check_as_datatype(MegaraDataType.IMAGE_FLAT)
    with pytest.raises(ValidationError, match="LAMPS3S"):
        checker(hdulist)
