#
# Copyright 2016-2025 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

"""Validators for Observing modes"""

import logging
import pkgutil
from io import StringIO

import json
import jsonschema

from numina.exceptions import ValidationError

from megaradrp.datatype import MegaraDataType

_logger = logging.getLogger(__name__)


def validate_focus(mode, obresult):
    """Validate FOCUS_SPECTROGRAPH"""
    image_groups = {}
    for idx, frame in enumerate(obresult.frames):
        with frame.open() as img:
            try:
                focus_val = img[0].header["FOCUS"]
                # FIXME: This should be done using an scheme for MEGARA
                if not isinstance(focus_val, (int, float)):
                    raise ValidationError(f"FOCUS must be integer, not {type(focus_val)}")
                if focus_val not in image_groups:
                    image_groups[focus_val] = []
                image_groups[focus_val].append(frame)
            except Exception as exc:
                raise ValidationError from exc

    if len(image_groups) < 2:
        raise ValidationError(f"We have only {len(image_groups)} different focus in OB")

    return True


def validate_key(mode, obresult, key):
    """Validate key"""

    # Assume that the individual images are valid IMG_COMP
    # check consistency of key
    kval = []
    for idx, frame in enumerate(obresult.frames):
        # SPECLAMP values
        with frame.open() as img:
            try:
                spec_val = img[0].header[key]
                kval.append(spec_val)
            except Exception as exc:
                raise ValidationError from exc

    if kval[:-1] == kval[1:]:
        return True
    else:
        raise ValidationError(f"{key} value is incorrect")


def validate_arc(mode, obresult):
    """Validate ARC_CALIBRATION"""

    # Assume that the individual images are valid IMG_COMP
    return validate_key(mode, obresult, "SPECLAMP")


def validate_flat(mode, obresult):
    """Validate FLAT"""

    # Assume that the individual images are valid IMAGE_FLAT
    return True


def validate_keyword_exists(header, key):
    """Verify that the keyword exists"""
    value = header.get(key)
    if value is None:
        msg = f'Expected keyword "{key}" is not present'
        raise ValidationError(msg)
    return True


def validate_keyword_value(header, key, expected):
    validate_keyword_exists(header, key)
    value = header.get(key)

    if value != expected:
        msg = f'Keyword "{key}" has value "{value}" != "{expected}"'
        raise ValidationError(msg)


def validate_keyword_any_value(header, key, any_expected):
    """Validate that keyword has any of allowed values"""
    validate_keyword_exists(header, key)
    value = header.get(key)

    for expected in any_expected:
        if value == expected:
            break
    else:
        msg = f'Keyword "{key}" has value "{value}" not in "{any_expected}"'
        raise ValidationError(msg)


def convert_headers(hdulist):
    headers = [convert_header(hdu.header) for hdu in hdulist]
    return headers


def convert_header(header):
    hdu_v = {}
    hdu_c = {}
    hdu_o = []
    hdu_repr = {"values": hdu_v, "comments": hdu_c, "ordering": hdu_o}

    for card in header.cards:
        key = card.keyword
        value = card.value
        comment = card.comment
        hdu_v[key] = value
        hdu_c[key] = comment
        hdu_o.append(key)

    return hdu_repr


def check_null(obj, level=None):
    return True


def check_invalid(obj, level=None):
    raise ValidationError


def _leaf_errors(error):
    if error.context:
        for sub_error in error.context:
            yield from _leaf_errors(sub_error)
    else:
        yield error


def relevant_error(error):
    """The most relevant leaf error of a schema error

    In an error of anyOf or oneOf, the relevant alternative is the
    one that fails deepest in the instance, with fewer errors.
    """
    if not error.context:
        return error
    alternatives = {}
    for sub_error in error.context:
        alternatives.setdefault(sub_error.relative_schema_path[0], []).extend(_leaf_errors(sub_error))
    leaves = min(alternatives.values(), key=lambda errors: (-max(len(e.absolute_path) for e in errors), len(errors)))
    return max(leaves, key=lambda e: len(e.absolute_path))


def validate_schema(validator, instance, hdu=None):
    """Validate `instance` with `validator`, raise ValidationError with a short message

    `hdu` is the index of the HDU when `instance` are the values of its header.
    """
    error = jsonschema.exceptions.best_match(validator.iter_errors(instance))
    if error is None:
        return
    error = relevant_error(error)
    path = list(error.absolute_path)
    if hdu is not None:
        path = [hdu, "values"] + path
    if len(path) >= 3 and path[1] == "values":
        where = f"HDU {path[0]}, keyword {path[2]}: "
    elif len(path) == 2 and path[1] == "values":
        where = f"HDU {path[0]}: "
    elif path:
        where = "/".join(str(part) for part in path) + ": "
    else:
        where = ""
    raise ValidationError(where + error.message)


class Checker:
    def __init__(self, validator):
        super(Checker, self).__init__()
        self.validator = validator

    def check(self, obj, level=None):

        return self.check_post(obj, level=level)

    def __call__(self, hdulist, level=None):
        return self.check(hdulist, level=level)

    def check_post(self, hdulist, level=None):
        return True


class ImageChecker(Checker):
    def __init__(self, validator):
        super(ImageChecker, self).__init__(validator)

    def check(self, hdulist, level=None):
        dheaders = convert_headers(hdulist)
        self.check_dheaders(dheaders, level=level)
        self.check_post(hdulist, level=level)
        return True

    def check_dheaders(self, dheaders, level=None):
        validate_schema(self.validator, dheaders)

        if len(dheaders) == 2:
            values_fibers = dheaders[1]["values"]
            values_primary = dheaders[0]["values"]
            check_header_additional(values_primary, values_fibers)


class StructChecker(Checker):
    def __init__(self, validator):
        super(StructChecker, self).__init__(validator)

    def check(self, obj, level=None):
        validate_schema(self.validator, obj)
        self.check_post(obj, level=level)
        return True


class BaseChecker(ImageChecker):
    def __init__(self, validator):
        super(BaseChecker, self).__init__(validator)


def definition_validator(validator, name):
    """Validator of the definition `name` of the schema of `validator`

    The other definitions are included, so that the references
    of the definition can be resolved.
    """
    schema = {"$ref": f"#/definitions/{name}", "definitions": validator.schema["definitions"]}
    return type(validator)(schema)


class ExtChecker(BaseChecker):
    def __init__(self, validator, definitions, n_ext=None):
        super(ExtChecker, self).__init__(validator)
        self.n_ext = n_ext
        self.definitions = definitions
        self.sub_validators = [definition_validator(validator, name) for name in definitions]

    def check_dheaders(self, dheaders, level=None):
        if self.n_ext is not None:
            if len(dheaders) != self.n_ext:
                msg = f"image has not expected number of HDUs ({self.n_ext})"
                raise ValueError(msg)

        super(ExtChecker, self).check_dheaders(dheaders, level=level)

        for sub_validator in self.sub_validators:
            validate_schema(sub_validator, dheaders[0]["values"], hdu=0)


def lamps_on(header, kind):
    """Keywords of the lamps of `kind` that are ON

    `kind` is "I" for the incandescent lamps and "S" for
    the comparison lamps.
    """
    nlamps = {"I": 2, "S": 5}[kind]
    keys = [f"LAMP{kind}{idx}S" for idx in range(1, nlamps + 1)]
    for key in keys:
        validate_keyword_exists(header, key)
    return [key for key in keys if header[key]]


class FlatImageChecker(ExtChecker):
    def __init__(self, validator):
        super(FlatImageChecker, self).__init__(validator, ["raw_hdu_values", "flat_hdu_values"], n_ext=2)

    def check_post(self, hdulist, level=None):
        """Check the lamps: some incandescent ON, all comparison OFF"""
        hdr = hdulist[0].header
        if not lamps_on(hdr, "I"):
            raise ValidationError("all incandescent lamps are OFF")
        if lamps_on(hdr, "S"):
            raise ValidationError("some comparison lamps are ON")


class CompImageChecker(ExtChecker):
    def __init__(self, validator):
        super(CompImageChecker, self).__init__(validator, ["raw_hdu_values", "comp_hdu_values"], n_ext=2)

    def check_post(self, hdulist, level=None):
        """Check the lamps: all incandescent OFF, some comparison ON"""
        hdr = hdulist[0].header
        if lamps_on(hdr, "I"):
            raise ValidationError("some incandescent lamps are ON")
        if not lamps_on(hdr, "S"):
            raise ValidationError("all comparison lamps are OFF")


class TargetImageChecker(ExtChecker):
    def __init__(self, validator):
        super(TargetImageChecker, self).__init__(validator, ["raw_hdu_values"], n_ext=2)

    def check_post(self, hdulist, level=None):
        """Check the lamps: all incandescent and comparison OFF"""
        hdr = hdulist[0].header
        if lamps_on(hdr, "I"):
            raise ValidationError("some incandescent lamps are ON")
        if lamps_on(hdr, "S"):
            raise ValidationError("some comparison lamps are ON")


class MasterFlatRSSChecker(ExtChecker):
    def __init__(self, validator):
        super(MasterFlatRSSChecker, self).__init__(validator, ["rss_naxis2"], n_ext=3)


class MasterSensitivityChecker(ExtChecker):
    def __init__(self, validator):
        super(MasterSensitivityChecker, self).__init__(validator, ["spec_hdu_values", "sensitivity_values"], n_ext=1)


def check_header_additional(values_primary, values_fibers):
    """Additional checks than can't be done with schema"""

    if values_primary["INSMODE"] != values_fibers["INSMODE"]:
        raise ValueError("insmode in PRIMARY != insmode in FIBERS")

    if values_fibers["INSMODE"] == "LCB":
        rbundles = [0, 93, 94, 95, 96, 97, 98, 99, 100]
    else:
        rbundles = range(1, 92 + 1)

    nfibers = values_fibers["NFIBERS"]
    # nbundles = values_fibers['NBUNDLES']

    for idbundle in rbundles:
        # types are check in the json schema
        for stype in ["P", "I", "T", "X", "Y", "O", "E"]:
            keyname = f"BUN{idbundle:03d}_{stype}"
            if keyname not in values_fibers:
                raise ValueError(f"keyname {keyname} not in values_fibers")

    missing_names = []
    for idfiber in range(1, nfibers + 1):
        # types are check in the json schema
        for stype in ["A", "D", "R", "X", "Y", "B"]:
            keyname = f"FIB{idfiber:03d}_{stype}"
            if keyname not in values_fibers:
                msg = f"keyname {keyname} not in values_fibers"
                raise ValueError(msg)

        # the names of the fibers are optional
        keyname = f"FIB{idfiber:03d}_N"
        if keyname not in values_fibers:
            missing_names.append(keyname)

    if missing_names:
        _logger.warning("%d keywords FIBxxx_N not in values_fibers", len(missing_names))


class CheckAsDatatype:
    """Collection of schemas for validation"""

    def __init__(self):

        image_schema_path = "baseimage.json"
        json_schema_path = "basestruct.json"

        data_image = pkgutil.get_data("megaradrp.schemas", image_schema_path)
        data_json = pkgutil.get_data("megaradrp.schemas", json_schema_path)
        schema_image = json.load(StringIO(data_image.decode("utf8")))
        schema_json = json.load(StringIO(data_json.decode("utf8")))

        ValClass = jsonschema.validators.validator_for(schema_image)
        self.validator_image = ValClass(schema_image)
        ValClass = jsonschema.validators.validator_for(schema_json)
        self.validator_json = ValClass(schema_json)

        raw_checker = ExtChecker(self.validator_image, ["raw_hdu_values"])
        proc_checker = ExtChecker(self.validator_image, ["proc_hdu_values"])
        rss_checker = ExtChecker(self.validator_image, ["rss_naxis2"])
        spec_checker = ExtChecker(self.validator_image, ["spec_hdu_values"])
        sens_checker = MasterSensitivityChecker(self.validator_image)
        struct_checker = StructChecker(self.validator_json)

        _megara_checkers = {}
        _megara_checkers[MegaraDataType.UNKNOWN] = check_null
        _megara_checkers[MegaraDataType.IMAGE_RAW] = raw_checker
        _megara_checkers[MegaraDataType.IMAGE_BIAS] = ExtChecker(
            self.validator_image, ["raw_hdu_values", "bias_hdu_values"], n_ext=1
        )
        _megara_checkers[MegaraDataType.IMAGE_DARK] = ExtChecker(
            self.validator_image, ["raw_hdu_values", "dark_hdu_values"], n_ext=1
        )
        _megara_checkers[MegaraDataType.IMAGE_SLITFLAT] = raw_checker
        _megara_checkers[MegaraDataType.IMAGE_FLAT] = FlatImageChecker(self.validator_image)
        _megara_checkers[MegaraDataType.IMAGE_COMP] = CompImageChecker(self.validator_image)
        #
        _megara_checkers[MegaraDataType.IMAGE_TWILIGHT] = raw_checker
        _megara_checkers[MegaraDataType.IMAGE_TEST] = raw_checker
        _megara_checkers[MegaraDataType.IMAGE_TARGET] = TargetImageChecker(self.validator_image)
        #
        _megara_checkers[MegaraDataType.IMAGE_PROCESSED] = proc_checker
        _megara_checkers[MegaraDataType.MASTER_BPM] = ExtChecker(
            self.validator_image, ["proc_hdu_values", "master_bpm_hdu_values"], n_ext=1
        )
        _megara_checkers[MegaraDataType.MASTER_BIAS] = ExtChecker(
            self.validator_image, ["proc_hdu_values", "master_bias_hdu_values"], n_ext=1
        )
        _megara_checkers[MegaraDataType.MASTER_DARK] = ExtChecker(
            self.validator_image, ["proc_hdu_values", "master_dark_hdu_values"], n_ext=1
        )
        _megara_checkers[MegaraDataType.MASTER_SLITFLAT] = proc_checker
        _megara_checkers[MegaraDataType.DIFFUSE_LIGHT] = proc_checker
        #
        _megara_checkers[MegaraDataType.RSS_PROCESSED] = rss_checker
        _megara_checkers[MegaraDataType.RSS_WL_PROCESSED] = rss_checker
        _megara_checkers[MegaraDataType.MASTER_FLAT] = MasterFlatRSSChecker(self.validator_image)
        _megara_checkers[MegaraDataType.MASTER_TWILIGHT] = rss_checker

        _megara_checkers[MegaraDataType.SPEC_PROCESSED] = spec_checker
        _megara_checkers[MegaraDataType.MASTER_SENSITIVITY] = sens_checker

        _megara_checkers[MegaraDataType.STRUCT_PROCESSED] = struct_checker
        _megara_checkers[MegaraDataType.TRACE_MAP] = struct_checker
        _megara_checkers[MegaraDataType.MODEL_MAP] = struct_checker
        _megara_checkers[MegaraDataType.WAVE_CALIB] = struct_checker

        self._megara_checkers = _megara_checkers

    def __call__(self, datatype):
        return self._megara_checkers[datatype]


check_as_datatype = CheckAsDatatype()
