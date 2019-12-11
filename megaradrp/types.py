#
# Copyright 2011-2019 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Products of the Megara Pipeline"""

import logging

from numina.core import DataFrameType
from numina.types.product import DataProductMixin
from numina.types.datatype import DataType
from numina.types.array import ArrayType
from numina.types.linescatalog import LinesCatalog
from numina.exceptions import ValidationError

from megaradrp.datamodel import MegaraDataModel, QueryAttribute


_logger = logging.getLogger(__name__)


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


def validate_fiber_ext(header_f):
    _logger.debug('validate fiber extension')


class MegaraFrame(DataFrameType):
    """A processed frame"""

    tags_headers = {}

    def __init__(self, *args, **kwds):
        super(MegaraFrame, self).__init__(datamodel=MegaraDataModel)

    def validate_hdulist(self, hdulist):
        # nhdus = len(hdulist)
        header_0 = hdulist[0].header

        validate_keyword_value(header_0, 'INSTRUME', 'MEGARA')

        kexits = ['DATE-OBS', 'UUID']
        # INSMODE can be LCB or MOS

        for key in kexits:
            validate_keyword_exists(header_0, key)

        tagsk = getattr(self, '__tags__', [])
        # FIXME: this is really ugly
        keymap = self.datamodel.extractor.map
        for key in tagsk:
            ksc = keymap[key]
            validate_keyword_exists(hdulist[ksc.ext].header, ksc.key)

    def validate_fibers(self, hdulist):
        if 'FIBERS' not in hdulist:
            msg = '"{}" extension not found'.format('FIBERS')
            raise ValidationError(hdulist, msg)

        header_f = hdulist['FIBERS'].header

        validate_fiber_ext(header_f)


class ProcessedFrame(MegaraFrame):
    """A processed frame"""

    tags_headers = {}


class ProcessedImage(ProcessedFrame):
    """A processed image"""

    def validate_hdulist(self, hdulist):

        super(ProcessedImage, self).validate_hdulist(hdulist)
        header_0 = hdulist[0].header
        validate_keyword_value(header_0, 'NAXIS1', 4096)
        validate_keyword_value(header_0, 'NAXIS2', 4112)


class ProcessedRSS(ProcessedFrame):
    """A processed RSS image"""

    def validate_hdulist(self, hdulist):
        super(ProcessedRSS, self).validate_hdulist(hdulist)
        _logger.debug('validate ProcessedRSS')

        header_0 = hdulist[0].header

        # validate_keyword_value(header_0, 'NUMTYPE', 'MasterBias')
        # 4096 for non WCS, and 4300 for WCS
        validate_keyword_any_value(header_0, 'NAXIS1', [4096, 4300])

        # This can be 623 or 644
        validate_keyword_any_value(header_0, 'NAXIS2', [623, 644])

        self.validate_fibers(hdulist)


class ProcessedMultiRSS(ProcessedFrame):
    """A processed RSS image not to be stored"""
    pass


class ProcessedSpectrum(ProcessedFrame):
    """A 1d spectrum"""
    pass


class ProcessedImageProduct(DataProductMixin, ProcessedImage):
    pass


class ProcessedRSSProduct(DataProductMixin, ProcessedRSS):
    pass


class ProcessedSpectrumProduct(DataProductMixin, ProcessedSpectrum):
    pass


class MegaraLinesCatalog(LinesCatalog):
    __tags__ = {'speclamp': QueryAttribute('speclamp', str), 'vph': QueryAttribute('vph', str)}
    # We are not passing the table of query_attrs in datamodel

    def name(self):
        return "LinesCatalog"


class MasterBias(ProcessedImageProduct):
    """A Master Bias image"""

    def validate_hdulist(self, hdulist):
        super(MasterBias, self).validate_hdulist(hdulist)
        _logger.debug('validate MasterBias')


class MasterTwilightFlat(ProcessedRSSProduct):
    __tags__ = ['insmode', 'vph', 'confid']

    def validate_hdulist(self, hdulist):
        super(MasterTwilightFlat, self).validate_hdulist(hdulist)
        _logger.debug('validate MasterTwilightFlat')


class MasterDark(ProcessedImageProduct):
    """A Master Dark image"""
    pass


class MasterFiberFlat(ProcessedRSSProduct):
    __tags__ = ['insmode', 'vph', 'confid']

    def validate_hdulist(self, hdulist):
        super(MasterFiberFlat, self).validate_hdulist(hdulist)
        _logger.debug('validate MasterFiberFlat')


class MasterSlitFlat(ProcessedImageProduct):
    __tags__ = ['insmode', 'vph']

    def validate_hdulist(self, hdulist):
        super(MasterSlitFlat, self).validate_hdulist(hdulist)
        _logger.debug('validate MasterSlitFlat')


class MasterBPM(ProcessedImageProduct):
    """Bad Pixel Mask product"""

    def validate_hdulist(self, hdulist):
        _logger.debug('validate MasterBPM')
        super(MasterBPM, self).validate_hdulist(hdulist)


class DiffuseLightCorrection(ProcessedImageProduct):
    """Image to correct from diffuse light"""
    pass


class MasterSensitivity(ProcessedSpectrumProduct):
    """Sensitivity correction."""
    pass


class ReferenceExtinctionTable(DataProductMixin, ArrayType):
    """Atmospheric Extinction."""

    def validate(self, obj):
        if obj is None:
            # None is valid
            pass
        else:
            super(ReferenceExtinctionTable, self).validate(obj)

    def convert(self, obj):
        # Support None value
        if obj is None:
            return None
        else:
            return super(ReferenceExtinctionTable, self).convert(obj)


class ReferenceSpectrumTable(DataProductMixin, ArrayType):
    """The spectrum of a reference star"""
    pass


class JSONstorage(DataType):
    def __init__(self, default=None):
        super(JSONstorage, self).__init__(ptype=dict, default=default)

    def _datatype_dump(self, obj, where):
        import json
        filename = where + '.json'

        with open(filename, 'w') as fd:
            fd.write(json.dumps(obj, sort_keys=True, indent=2,
                                separators=(',', ': ')))

        return filename

    def _datatype_load(self, obj):
        import json
        try:
            with open(obj, 'r') as fd:
                data = json.load(fd)
        except IOError as e:
            raise e
        return data


class FocusWavelength(JSONstorage):
    """Rich table with focus and wavelength"""
    pass
