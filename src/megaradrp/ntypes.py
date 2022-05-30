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
from numina.types.array import ArrayType, ArrayNType
from numina.types.linescatalog import LinesCatalog
from numina.exceptions import ValidationError

import megaradrp.validators as valid
from megaradrp.datatype import MegaraDataType
from megaradrp.datamodel import MegaraDataModel, QueryAttribute


_logger = logging.getLogger(__name__)


def validate_fiber_ext(header_f):
    _logger.debug('validate fiber extension')


class Point2D(ArrayNType):
    """A type of 2D cartesian coordinates"""
    def __init__(self, *args, **kwds):
        super(ArrayNType, self).__init__(2)


class MegaraFrame(DataFrameType):
    """A processed frame"""
    DATATYPE = MegaraDataType.IMAGE_RAW
    tags_headers = {}

    def __init__(self, *args, **kwds):
        super(MegaraFrame, self).__init__(datamodel=MegaraDataModel)

    def validate_hdulist(self, hdulist):
        _logger.debug('validate MasterBias')
        checker = valid.check_as_datatype(self.DATATYPE)
        return checker(hdulist)

    def extract_tags(self, obj):
        """Extract tags from serialized file"""

        objl = self.convert(obj)
        ext = self.datamodel.extractor_map['fits']
        tags = {}

        if objl:
            with objl.open() as hdulist:
                for field in self.names_t:
                    tags[field] = ext.extract(field, hdulist)
                return tags
        else:
            return {}


class ProcessedFrame(MegaraFrame):
    """A processed frame"""
    DATATYPE = MegaraDataType.IMAGE_PROCESSED
    tags_headers = {}


class ProcessedImage(ProcessedFrame):
    """A processed image"""
    DATATYPE = MegaraDataType.IMAGE_PROCESSED


class ProcessedRSS(ProcessedFrame):
    """A processed RSS image"""
    DATATYPE = MegaraDataType.RSS_PROCESSED


class ProcessedMultiRSS(ProcessedFrame):
    """A processed RSS image not to be stored"""
    pass


class ProcessedSpectrum(ProcessedFrame):
    """A 1d spectrum"""
    DATATYPE = MegaraDataType.SPEC_PROCESSED
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
    DATATYPE = MegaraDataType.MASTER_BIAS


class MasterTwilightFlat(ProcessedRSSProduct):
    __tags__ = ['insmode', 'vph', 'confid']
    DATATYPE = MegaraDataType.MASTER_TWILIGHT


class MasterDark(ProcessedImageProduct):
    """A Master Dark image"""
    DATATYPE = MegaraDataType.MASTER_DARK


class MasterFiberFlat(ProcessedRSSProduct):
    __tags__ = ['insmode', 'vph', 'confid']
    DATATYPE = MegaraDataType.MASTER_FLAT


class MasterSlitFlat(ProcessedImageProduct):
    __tags__ = ['insmode', 'vph']
    DATATYPE = MegaraDataType.MASTER_SLITFLAT


class MasterBPM(ProcessedImageProduct):
    """Bad Pixel Mask product"""
    DATATYPE = MegaraDataType.MASTER_BPM


class DiffuseLightCorrection(ProcessedImageProduct):
    """Image to correct from diffuse light"""
    pass


class MasterSensitivity(ProcessedSpectrumProduct):
    """Sensitivity correction."""
    pass


class SkyRSS(ProcessedRSS):
    """A processed RSS image"""
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
