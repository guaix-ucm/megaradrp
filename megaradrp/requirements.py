#
# Copyright 2015-2019 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Typical requirements of recipes"""

from numina.core import Requirement
from numina.types.multitype import MultiType

import megaradrp.ntypes
import megaradrp.products
import megaradrp.products.modelmap


class MasterBiasRequirement(Requirement):
    """Requirement for Master Bias images"""
    def __init__(self, optional=False):
        super(MasterBiasRequirement, self).__init__(
            megaradrp.ntypes.MasterBias,
            'Master BIAS image',
            optional=optional,
            validation=True
        )


class MasterBPMRequirement(Requirement):
    def __init__(self, optional=True):
        super(MasterBPMRequirement, self).__init__(
            megaradrp.ntypes.MasterBPM,
            'Master Bad Pixel Mask',
            optional=optional,
            validation=False
        )


class MasterDarkRequirement(Requirement):
    def __init__(self, optional=True):
        super(MasterDarkRequirement, self).__init__(
            megaradrp.ntypes.MasterDark,
            'Master DARK image',
            optional=optional,
            validation=True
        )


class MasterFiberFlatRequirement(Requirement):
    def __init__(self):
        super(MasterFiberFlatRequirement, self).__init__(
            megaradrp.ntypes.MasterFiberFlat,
            'Master fiber flat calibration',
            validation=True
        )


class MasterSlitFlatRequirement(Requirement):
    def __init__(self, optional=True):
        super(MasterSlitFlatRequirement,
              self).__init__(megaradrp.ntypes.MasterSlitFlat,
                             'Master slit flat calibration',
                             optional=optional
                             )


class MasterTwilightRequirement(Requirement):
    def __init__(self, optional=True):
        super(MasterTwilightRequirement,
              self).__init__(megaradrp.ntypes.MasterTwilightFlat,
                             'Master twlight flat calibration',
                             optional=optional
                             )


class MasterTraceMapRequirement(Requirement):
    def __init__(self):
        super(MasterTraceMapRequirement, self).__init__(
            megaradrp.products.TraceMap,
            'Trace information of the Apertures',
            validation=True
        )


class MasterAperturesRequirement(Requirement):
    def __init__(self, alias=None):
        super(MasterAperturesRequirement, self).__init__(MultiType(
            megaradrp.products.modelmap.ModelMap,
            megaradrp.products.TraceMap),
            'Apertures information for extraction',
            validation=True,
            alias=alias
        )


class WavelengthCalibrationRequirement(Requirement):
    def __init__(self):
        super(WavelengthCalibrationRequirement, self).__init__(
            megaradrp.products.WavelengthCalibration,
            'Wavelength calibration table',
            validation=True
        )


class LinesCatalogRequirement(Requirement):
    def __init__(self):
        super(LinesCatalogRequirement, self).__init__(megaradrp.ntypes.MegaraLinesCatalog, 'Catalog of lines')


class SkyRSSRequirement(Requirement):
    def __init__(self, optional=True):
        super(SkyRSSRequirement, self).__init__(
            megaradrp.ntypes.SkyRSS,
            'Row Stacked Spectra of the sky',
            optional=optional
        )


class SensitivityRequirement(Requirement):
    def __init__(self, optional=True):
        super(SensitivityRequirement,
              self).__init__(megaradrp.ntypes.MasterSensitivity,
                             'Master sensitivity for flux calibration',
                             optional=optional
                             )


class ReferenceExtinction(Requirement):
    def __init__(self, optional=True):
        super(ReferenceExtinction,
              self).__init__(megaradrp.ntypes.ReferenceExtinctionTable,
                             "Reference extinction",
                             optional=optional
                             )


class DiffuseLightRequirement(Requirement):
    def __init__(self, optional=True):
        super(DiffuseLightRequirement,
              self).__init__(megaradrp.ntypes.DiffuseLightCorrection,
                             'Diffuse light correction image',
                             optional=optional
                             )
