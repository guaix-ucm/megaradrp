#
# Copyright 2016-2020 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Validators for Observing modes"""

import sys
import six
import pkgutil
from six import StringIO

import json
import jsonschema

from numina.exceptions import ValidationError

from megaradrp.datatype import MegaraDataType


def validate_focus(mode, obresult):
    """Validate FOCUS_SPECTROGRAPH"""
    image_groups = {}
    for idx, frame in enumerate(obresult.frames):
        with frame.open() as img:
            try:
                focus_val = img[0].header['FOCUS']
                # FIXME: This should be done using an scheme for MEGARA
                if not isinstance(focus_val, (int, float)):
                    raise ValidationError("FOCUS must be integer, not {}".format(type(focus_val)))
                if focus_val not in image_groups:
                    image_groups[focus_val] = []
                image_groups[focus_val].append(frame)
            except Exception:
                _type, exc, tb = sys.exc_info()
                six.reraise(ValidationError, exc, tb)

    if len(image_groups) < 2:
        raise ValidationError('We have only {} different focus in OB'.format(len(image_groups)))

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
            except Exception:
                _type, exc, tb = sys.exc_info()
                six.reraise(ValidationError, exc, tb)

    if kval[:-1] == kval[1:]:
        return True
    else:
        raise ValidationError("{} value is incorrect".format(key))


def validate_arc(mode, obresult):
    """Validate ARC_CALIBRATION"""

    # Assume that the individual images are valid IMG_COMP
    return validate_key(mode, obresult, 'SPECLAMP')


def validate_flat(mode, obresult):
    """Validate FLAT"""

    # Assume that the individual images are valid IMG_COMP
    return True


def validate_keyword_exists(header, key):
    """Verify that the keyword exists"""
    value = header.get(key)
    if value is None:
        msg = 'Expected keyword "{}" is not present'.format(key)
        raise ValidationError(msg)
    return True


def validate_keyword_value(header, key, expected):
    from numina.exceptions import ValidationError

    validate_keyword_exists(header, key)
    value = header.get(key)

    if value != expected:
        msg = 'Keyword "{0}" has value "{1}" != "{2}"'.format(key, value, expected)
        raise ValidationError(msg)


def validate_keyword_any_value(header, key, any_expected):
    """Validate that keyword has any of allowed values"""
    from numina.exceptions import ValidationError

    validate_keyword_exists(header, key)
    value = header.get(key)

    for expected in any_expected:
        if value == expected:
            break
    else:
        msg = 'Keyword "{0}" has value "{1}" not in "{2}"'.format(key, value, any_expected)
        raise ValidationError(msg)


def convert_headers(hdulist):
    headers = [convert_header(hdu.header) for hdu in hdulist]
    return headers


def convert_header(header):
    hdu_v = {}
    hdu_c = {}
    hdu_o = []
    hdu_repr = {'values': hdu_v, 'comments': hdu_c, 'ordering': hdu_o}

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


class Checker(object):
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
        super(ImageChecker, self).check_post(hdulist, level=level)
        return True

    def check_dheaders(self, dheaders, level=None):
        # Check with json schema
        # convert header to dict
        try:
            self.validator.validate(dheaders)
        except jsonschema.exceptions.ValidationError:
            raise

        if len(dheaders) == 2:
            values_fibers = dheaders[1]['values']
            values_primary = dheaders[0]['values']
            check_header_additional(values_primary, values_fibers)


class StructChecker(Checker):
    def __init__(self, validator):
        super(StructChecker, self).__init__(validator)

    def check(self, obj, level=None):
        self.validator.validate(obj)
        self.check_post(obj, level=level)
        return True


class BaseChecker(ImageChecker):
    def __init__(self, validator):
        super(BaseChecker, self).__init__(validator)


# TODO: insert all subschemas in the general schema
_sub_schema_rss = {
    "oneOf": [
    {
        "type": "object",
        "properties": {
            #"NAXIS1": {"const": 4300},
            "NAXIS2": {"const": 623},
            "INSMODE": {"const": "LCB"}
        }
    },
    {
        "type": "object",
        "properties": {
            #"NAXIS1": {"const": 4300},
            "NAXIS2": {"const": 644},
            "INSMODE": {"const": "MOS"}
        }
    }
    ]
}

_sub_schema_bias = {
    "type": "object",
    "properties": {
        "OBJECT": {"const": "BIAS"},
        "OBSMODE": {"const": "MegaraBiasImage"},
        "IMAGETYP": {"const": "IMAGE_BIAS"},
        "EXPTIME": {"type": "number", "maximum": 0},
        "DARKTIME": {"type": "number", "maximum": 0}
    }
}


_sub_schema_master_bpm = {
    "type": "object",
    "properties": {
        "NUMTYPE": {"enum": ['MasterBias', 'MASTER_BPM']}
    }
}

_sub_schema_master_bias = {
    "type": "object",
    "properties": {
        #"OBJECT": {"const": "BIAS"},
        "OBSMODE": {"const": "MegaraBiasImage"},
        "IMAGETYP": {"const": "MASTER_BIAS"},
        "EXPTIME": {"type": "number", "maximum": 0},
        "DARKTIME": {"type": "number", "maximum": 0},
        "NUMTYPE": {"enum": ['MasterBias', 'MASTER_BIAS']}
    }
}


_sub_schema_dark = {
    "type": "object",
    "properties": {
        "OBSMODE": {"const": "MegaraDarkImage"},
        "IMAGETYP": {"const": "IMAGE_DARK"},
    }
}


_sub_schema_master_dark = {
    "type": "object",
    "properties": {
        "OBSMODE": {"const": "MegaraDarkImage"},
        "IMAGETYP": {"const": "MASTER_DARK"},
    }
}


class ExtChecker(BaseChecker):
    def __init__(self, schema, sub_schemas, n_ext=None):
        super(ExtChecker, self).__init__(schema)
        self.n_ext = n_ext
        self.sub_schemas = sub_schemas

    def check_dheaders(self, dheaders, level=None):

        try:
            super(ExtChecker, self).check_dheaders(dheaders, level=level)
        except jsonschema.exceptions.ValidationError:
            pass
        # Image must have only one extension
        if self.n_ext is not None:
            if len(dheaders) != self.n_ext:
                msg = 'image has not expected number of HDUs ({})'.format(self.n_ext)
                raise ValueError(msg)

        for sub_schema in self.sub_schemas:
            try:
                if isinstance(sub_schema, str):
                    url, fragment = self.validator.resolver.resolve(sub_schema)
                else:
                    fragment = sub_schema

                jsonschema.validate(dheaders[0]['values'], schema=fragment)
            except jsonschema.exceptions.ValidationError:
                raise


class FlatImageChecker(ExtChecker):
    def __init__(self, schema):

        _sub_schema_flat = {
            "type": "object",
            "properties": {
                "OBSMODE": {"enum": [
                    "MegaraFiberFlatImage", "MegaraTraceMap", "MegaraModelMap", "MegaraSuccess"]
                },
                "IMAGETYP": {"const": "IMAGE_FLAT"},
            }
        }

        super(FlatImageChecker, self).__init__(schema, ["#/definitions/raw_hdu_values", _sub_schema_flat], n_ext=2)

    def check_post(self, hdulist, level=None):
        """Additional checks"""
        self.check_post_level1(hdulist)

    def check_post_level1(self, hdulist):
        """Additional checks"""
        hdr = hdulist[0].header
        # Flat must have inc LAMPS-ON
        # Flat must have comp LAMPS-OFF
        lamp_i_s = (hdr['LAMPI1S'] or hdr['LAMPI1S'])
        if not lamp_i_s:
            msg = 'all incandescent lamps are OFF'
            raise ValidationError(msg)
        lamp_s_s = True
        for idx in range(1, 6):
            label = 'LAMPS{}S'.format(idx)
            lamp_s_s = lamp_s_s and hdr[label]
        if lamp_s_s:
            msg = 'some comparation lamps are ON'
            raise ValidationError(msg)


class CompImageChecker(ExtChecker):
    def __init__(self, schema):

        _sub_schema_comp = {
            "type": "object",
            "properties": {
                "OBSMODE": {"enum": ["MegaraArcCalibration", "MegaraSuccess"]},
                "IMAGETYP": {"const": "IMAGE_COMP"},
            }
        }

        super(CompImageChecker, self).__init__(schema, ["#/definitions/raw_hdu_values", _sub_schema_comp], n_ext=2)

    def check_post(self, hdulist, level=None):
        """Additional checks"""
        self.check_post_level1(hdulist)

    def check_post_level1(self, hdulist):
        """Additional checks"""
        hdr = hdulist[0].header
        # Flat must have all inc LAMPS-OFF
        # Flat must have some comp LAMPS-ON
        lamp_i_s = (hdr['LAMPI1S'] or hdr['LAMPI1S'])
        if lamp_i_s:
            msg = 'some incandescent lamps are ON'
            raise ValidationError(msg)
        lamp_s_s = False
        for idx in range(1, 6):
            label = 'LAMPS{}S'.format(idx)
            lamp_s_s = lamp_s_s or hdr[label]
        if not lamp_s_s:
            msg = 'all comparation lamps are OFF'
            raise ValidationError(msg)


class TargetImageChecker(ExtChecker):
    def __init__(self, schema):

        _sub_schema_target = {
            "type": "object",
            "properties": {
            }
        }

        super(TargetImageChecker, self).__init__(schema,
            ["#/definitions/raw_hdu_values", _sub_schema_target], n_ext=2
        )

    def check_post(self, hdulist, level=None):
        """Additional checks"""
        self.check_post_level1(hdulist)

    def check_post_level1(self, hdulist):
        """Additional checks"""
        hdr = hdulist[0].header
        # Target must have all inc LAMPS-OFF
        # Target must have all comp LAMPS-ON
        lamp_i_s = (hdr['LAMPI1S'] or hdr['LAMPI1S'])
        if lamp_i_s:
            msg = 'some incandescent lamps are ON'
            raise ValidationError(msg)
        lamp_s_s = False
        for idx in range(1, 6):
            label = 'LAMPS{}S'.format(idx)
            lamp_s_s = lamp_s_s or hdr[label]
        if lamp_s_s:
            msg = 'some comparation lamps are ON'
            raise ValidationError(msg)


class MasterFlatRSSChecker(ExtChecker):
    def __init__(self, schema):
        super(MasterFlatRSSChecker, self).__init__(schema,
            [_sub_schema_rss], n_ext=3
        )

class MasterSensitivityChecker(ExtChecker):
    def __init__(self, schema):
        super(MasterSensitivityChecker, self).__init__(schema,
            ["#/definitions/spec_hdu_values", "#/definitions/sensitivity_values"], n_ext=1
        )


def check_header_additional(values_primary, values_fibers):
    """Additional checks than can't be done with schema"""


    if values_primary['INSMODE'] != values_fibers['INSMODE']:
        raise ValueError('insmode in PRIMARY != insmode in FIBERS')

    if values_fibers['INSMODE'] == 'LCB':
        rbundles = [0, 93, 94, 95, 96, 97, 98, 99, 100]
    else:
        rbundles = range(1, 92 + 1)

    nfibers = values_fibers['NFIBERS']
    nbundles = values_fibers['NBUNDLES']

    for idbundle in rbundles:
        # types are check in the json schema
        for stype in ['P', "I", "T", "X", "Y", "O", "E"]:
            keyname = "BUN{:03d}_{}".format(idbundle, stype)
            if keyname not in values_fibers:
                raise ValueError("keyname {} not in values_fibers".format(keyname))

    for idfiber in range(1, nfibers + 1):
        # types are check in the json schema
        for stype in ['A', "D", "R", "X", "Y", "B"]:
            keyname = "FIB{:03d}_{}".format(idfiber, stype)
            if keyname not in values_fibers:
                msg = "keyname {} not in values_fibers".format(keyname)
                raise ValueError(msg)

        for stype in ["N"]:
            keyname = "FIB{:03d}_{}".format(idfiber, stype)
            if keyname not in values_fibers:
                msg = "keyword {} not in values_fibers".format(keyname)
                print(msg)
                #raise ValueError(msg)


class CheckAsDatatype(object):
    """Collection of schemas for validation"""
    def __init__(self):

        image_schema_path = "baseimage.json"
        json_schema_path = "basestruct.json"

        data_image = pkgutil.get_data('megaradrp.schemas', image_schema_path)
        data_json = pkgutil.get_data('megaradrp.schemas', json_schema_path)
        schema_image = json.load(StringIO(data_image.decode('utf8')))
        schema_json = json.load(StringIO(data_json.decode('utf8')))

        ValClass = jsonschema.validators.validator_for(schema_image)
        self.validator_image = ValClass(schema_image)
        ValClass = jsonschema.validators.validator_for(schema_json)
        self.validator_json = ValClass(schema_json)

        raw_checker = ExtChecker(self.validator_image, ["#/definitions/raw_hdu_values"])
        proc_checker = ExtChecker(self.validator_image, ["#/definitions/proc_hdu_values"])
        rss_checker = ExtChecker(self.validator_image, [_sub_schema_rss])
        spec_checker = ExtChecker(self.validator_image, ["#/definitions/spec_hdu_values"])
        sens_checker = MasterSensitivityChecker(self.validator_image)
        struct_checker = StructChecker(self.validator_json)

        _megara_checkers = {}
        _megara_checkers[MegaraDataType.UNKNOWN] = check_null
        _megara_checkers[MegaraDataType.IMAGE_RAW] = raw_checker
        _megara_checkers[MegaraDataType.IMAGE_BIAS] = ExtChecker(self.validator_image, ["#/definitions/raw_hdu_values", _sub_schema_bias], n_ext=1)
        _megara_checkers[MegaraDataType.IMAGE_DARK] = ExtChecker(self.validator_image, ["#/definitions/raw_hdu_values", _sub_schema_dark], n_ext=1)
        _megara_checkers[MegaraDataType.IMAGE_SLITFLAT] = raw_checker
        _megara_checkers[MegaraDataType.IMAGE_FLAT] = FlatImageChecker(self.validator_image)
        _megara_checkers[MegaraDataType.IMAGE_COMP] = CompImageChecker(self.validator_image)
        #
        _megara_checkers[MegaraDataType.IMAGE_TWILIGHT] = raw_checker
        _megara_checkers[MegaraDataType.IMAGE_TEST] = raw_checker
        _megara_checkers[MegaraDataType.IMAGE_TARGET] = TargetImageChecker(self.validator_image)
        #
        _megara_checkers[MegaraDataType.IMAGE_PROCESSED] = proc_checker
        _megara_checkers[MegaraDataType.MASTER_BPM] = ExtChecker(self.validator_image, [
            "#/definitions/proc_hdu_values", _sub_schema_master_bpm
        ], n_ext=1)
        _megara_checkers[MegaraDataType.MASTER_BIAS] = ExtChecker(self.validator_image, ["#/definitions/proc_hdu_values", _sub_schema_master_bias], n_ext=1)
        _megara_checkers[MegaraDataType.MASTER_DARK] = ExtChecker(self.validator_image, ["#/definitions/proc_hdu_values", _sub_schema_master_dark], n_ext=1)
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
