#
# Copyright 2011-2017 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Products of the Megara Pipeline"""


from numina.core import DataFrameType
from numina.types.product import DataProductMixin
from numina.types.datatype import DataType
from numina.types.array import ArrayType
from numina.types.linescatalog import LinesCatalog


from megaradrp.datamodel import MegaraDataModel, QueryAttribute


class MegaraFrame(DataFrameType):
    """A processed frame"""

    tags_headers = {}

    def __init__(self, *args, **kwds):
        super(MegaraFrame, self).__init__(datamodel=MegaraDataModel)


class ProcessedFrame(MegaraFrame):
    """A processed frame"""

    tags_headers = {}


class ProcessedImage(ProcessedFrame):
    """A processed image"""
    pass


class ProcessedRSS(ProcessedFrame):
    """A processed RSS image"""
    pass


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
    pass


class MasterTwilightFlat(ProcessedRSSProduct):
    __tags__ = ['insmode', 'vph', 'confid']


class MasterDark(ProcessedImageProduct):
    """A Master Dark image"""
    pass


class MasterFiberFlat(ProcessedRSSProduct):
    __tags__ = ['insmode', 'vph', 'confid']


class MasterSlitFlat(ProcessedImageProduct):
    __tags__ = ['insmode', 'vph']


class MasterBPM(ProcessedImageProduct):
    """Bad Pixel Mask product"""
    pass


class MasterSensitivity(ProcessedSpectrumProduct):
    """Sensitivity correction."""
    pass


class ReferenceExtinctionTable(DataProductMixin, ArrayType):
    """Atmospheric Extinction."""
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
