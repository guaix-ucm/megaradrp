#
# Copyright 2016-2021 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Data model for MEGARA"""

from __future__ import division


import pkgutil
import logging

import astropy.io.fits as fits
from io import StringIO
from numina.datamodel import DataModel, QueryAttribute, KeyDefinition
from numina.util.convert import convert_date

from megaradrp.datatype import MegaraDataType, DataOrigin
import megaradrp.instrument.constants as cons


_logger = logging.getLogger(__name__)


class MegaraDataModel(DataModel):
    """Data model of MEGARA images"""

    query_attrs = {
        'vph': QueryAttribute('vph', str),
        'insmode': QueryAttribute('insmode', str),
        'insconf': QueryAttribute('insconf', str),
        'speclamp': QueryAttribute('speclamp', str),
        'temp': QueryAttribute('temp', float),
        'confid': QueryAttribute('confid', str)
    }

    meta_info_headers = [
        'instrument',
        'object',
        'observation_date',
        'uuid',
        'type',
        'mode',
        'exptime',
        'darktime',
        'insconf',
        'blckuuid',
        'quality_control',
        'vph',
        'insmode'
    ]

    db_info_keys = [
        'instrument',
        'object',
        'observation_date',
        'uuid',
        'type',
        'mode',
        'exptime',
        'darktime',
        'insconf',
        'blckuuid',
        'quality_control',
        'vph',
        'insmode'
    ]

    db_info_keys_extra = [
        'vph',
        'insmode'
    ]

    meta_dinfo_headers = [
        'exptime',
        'observation_date',
        'vph',
        'vphpos',
        'insmode',
        'focus',
        'osfilter',
        'uuid',
        'temp',
        'block_uuid',
        'insconf_uuid',
        'speclamp',
        'imgid'
    ]

    PLATESCALE = cons.GTC_FC_A_PLATESCALE.value

    def __init__(self):

        instrument_mappings = {
            'date_obs': ('DATE-OBS', 0, convert_date),
            'insconf': 'insconf',
            'insconf_uuid': 'insconf',
            'blckuuid': 'blckuuid',
            'block_uuid': 'blckuuid',
            'vph': ('VPH', 'undefined'),
            'vphpos': ('VPHWHPOS', 'undefined'),
            'focus': ('FOCUS', 'undefined'),
            'osfilter': ('OSFILTER', 'undefined'),
            'temp': ('SENTEMP4', 0.0),
            'speclamp': ('SPECLAMP', 'undefined'),
            'confid': KeyDefinition('CONFID', ext='FIBERS'),
        }
        super(MegaraDataModel, self).__init__(
            'MEGARA',
            instrument_mappings
        )

    def get_fiberconf(self, img):
        """Obtain FiberConf from image"""
        return get_fiberconf(img)

    def gather_info_oresult(self, val):
        return [self.gather_info_dframe(f) for f in val.images]

    def fiber_scale_unit(self, img, unit=False):
        return fiber_scale_unit(img, unit=unit)


def fiber_scale_unit(img, unit=False):
    funit = img['FIBERS'].header.get("FUNIT", "arcsec")

    if funit == "arcsec":
        scale = 1
    else:
        scale = cons.GTC_FC_A_PLATESCALE.value
    if unit:
        return scale, funit
    else:
        return scale


def get_fiberconf(img):
    """Obtain FiberConf from image"""

    main_insmode = img[0].header.get('INSMODE', '')

    if 'FIBERS' in img:
        # We have a 'fibers' extension
        # Information os there
        hdr_fiber = img['FIBERS'].header
        return read_fibers_extension(hdr_fiber, insmode=main_insmode)
    else:
        return get_fiberconf_default(main_insmode)


def create_default_fiber_header(insmode):
    """Obtain default FIBER header"""
    if insmode == 'LCB':
        slit_file = 'lcb_default_header.txt'
    elif insmode == 'MOS':
        slit_file = 'mos_default_header.txt'
    else:
        # Read fiber info from headers
        raise ValueError(f'Invalid INSMODE {insmode}')

    data = pkgutil.get_data('megaradrp.instrument.configs', slit_file)
    default_hdr = StringIO(data.decode('utf8'))
    hdr_fiber = fits.header.Header.fromfile(default_hdr)
    return hdr_fiber


def get_fiberconf_default(insmode):
    """Obtain default FiberConf object"""
    hdr_fiber = create_default_fiber_header(insmode)
    return read_fibers_extension(hdr_fiber)


def read_fibers_extension(hdr, insmode='LCB'):
    """Read the FIBERS extension

    Parameters
    ==========
    hdr :
       FITS header
    insmode : str
        default INSMODE

    Returns
    =======
    FibersConf


    """
    import megaradrp.instrument.focalplane as fp
    return fp.FocalPlaneConf.from_header(hdr)


def describe_hdulist_megara(hdulist):
    prim = hdulist[0].header
    instrument = prim.get("INSTRUME", "unknown")
    image_type = prim.get("IMAGETYP")
    if image_type is None:
        # try this also
        image_type = prim.get("NUMTYPE")

    # date_obs = convert_date(prim.get("DATE-OBS"))
    date_obs = prim.get("DATE-OBS")
    img_uuid = prim.get("UUID")
    insconf = prim.get("INSCONF", 'undefined')

    if image_type is None:
        # inferr from header
        datatype = megara_inferr_datetype_from_image(hdulist)
    else:
        datatype = MegaraDataType[image_type]

    obs = {}
    proc = {}
    if datatype.value < MegaraDataType.IMAGE_PROCESSED.value:
        origin = DataOrigin.OBSERVED
    else:
        origin = DataOrigin.PROCESSED

    return {'instrument': instrument, 'datatype': datatype,
            'origin': origin, 'uuid': img_uuid,
            'insconf': insconf,
            'observation': obs,
            'observation_date': date_obs,
            'processing': proc
            }


def megara_inferr_datatype(obj):

    if isinstance(obj, fits.HDUList):
        return megara_inferr_datetype_from_image(obj)
    elif isinstance(obj, dict):
        return megara_inferr_datetype_from_dict(obj)
    else:
        raise TypeError(f"I don't know how to inferr datatype from {obj}")


def megara_inferr_datetype_from_dict(obj):
    # this comes from JSON
    dtype = obj['type_fqn']
    if dtype in ["megaradrp.products.tracemap.TraceMap"]:
        return MegaraDataType.TRACE_MAP
    elif dtype in ["megaradrp.products.modelmap.ModelMap"]:
        return MegaraDataType.MODEL_MAP
    elif dtype in ["megaradrp.products.wavecalibration.WavelengthCalibration"]:
        return MegaraDataType.WAVE_CALIB
    else:
        return MegaraDataType.UNKNOWN


def megara_inferr_datetype_from_image(hdulist):
    IMAGE_RAW_SHAPE = (4212, 4196)
    IMAGE_PROC_SHAPE = (4112, 4096)
    RSS_IFU_PROC_SHAPE = (623, 4096)
    RSS_MOS_PROC_SHAPE = (644, 4096)
    RSS_IFU_PROC_WL_SHAPE = (623, 4300)
    RSS_MOS_PROC_WL_SHAPE = (644, 4300)
    SPECTRUM_PROC_SHAPE = (4300,)
    prim = hdulist[0].header

    image_type = prim.get("IMAGETYP")
    if image_type is None:
        # try this also
        image_type = prim.get("NUMTYPE")

    if image_type is not None:
        datatype = MegaraDataType[image_type]
        return datatype

    pshape = hdulist[0].shape
    obsmode = prim.get("OBSMODE", "unknown")
    if pshape == IMAGE_RAW_SHAPE:
        datatype = MegaraDataType.IMAGE_RAW
    elif pshape == IMAGE_PROC_SHAPE:
        datatype = MegaraDataType.IMAGE_PROCESSED
    elif pshape == RSS_IFU_PROC_SHAPE: # IFU
        datatype = MegaraDataType.RSS_PROCESSED
    elif pshape == RSS_MOS_PROC_SHAPE: # MOS
        datatype = MegaraDataType.RSS_PROCESSED
    elif pshape == RSS_IFU_PROC_WL_SHAPE: # IFU
        datatype = MegaraDataType.RSS_WL_PROCESSED
    elif pshape == RSS_MOS_PROC_WL_SHAPE: # MOS
        datatype = MegaraDataType.RSS_WL_PROCESSED
    elif pshape == SPECTRUM_PROC_SHAPE:
        datatype = MegaraDataType.SPEC_PROCESSED
    else:
        datatype = MegaraDataType.UNKNOWN

    if datatype == MegaraDataType.IMAGE_RAW:
        if obsmode in ["MegaraSuccess", "MegaraFail"]:
            sub_datatype = MegaraDataType.IMAGE_TEST
        elif obsmode in ["MegaraBiasImage"]:
            sub_datatype = MegaraDataType.IMAGE_BIAS
        elif obsmode in ["MegaraDarkImage"]:
            sub_datatype = MegaraDataType.IMAGE_DARK
        elif obsmode in ["MegaraSlitFlat"]:
            sub_datatype = MegaraDataType.IMAGE_SLITFLAT
        elif obsmode in ["MegaraFiberFlatImage", "MegaraTraceMap", "MegaraModelMap"]:
            sub_datatype = MegaraDataType.IMAGE_FLAT
        elif obsmode in ["MegaraArcCalibration"]:
            sub_datatype = MegaraDataType.IMAGE_COMP
        elif obsmode in ["MegaraTwilightFlatImage"]:
            sub_datatype = MegaraDataType.IMAGE_TWILIGHT
        elif obsmode in ["MegaraLcbImage", "MegaraMosImage"]:
            sub_datatype = MegaraDataType.IMAGE_TARGET
        elif obsmode in ["MegaraMosStdStar", "MegaraExtinctionStar",
                         "MegaraLcbStdStar", "MegaraSensitivityStar"]:
            sub_datatype = MegaraDataType.IMAGE_TARGET
        elif obsmode in ["MegaraFocusTelescope",
                         "MegaraLcbAcquisition", "MegaraMosAcquisition"]:
            sub_datatype = MegaraDataType.IMAGE_TARGET
        elif obsmode in ["MegaraBadPixelMask", "MegaraFocusSpectrograph"]:
            sub_datatype = MegaraDataType.IMAGE_RAW
        else:
            sub_datatype = MegaraDataType.UNKNOWN

        return sub_datatype
    elif datatype == MegaraDataType.SPEC_PROCESSED:
        numrnam = prim.get("NUMRNAM", "unknown")
        if numrnam in ['LCBStandardRecipe', 'MOSStandardRecipe']:
            sub_datatype = MegaraDataType.MASTER_SENSITIVITY
        else:
            sub_datatype = datatype
        return sub_datatype
    return datatype


def check_obj_megara(obj, astype=None, level=None):
    import megaradrp.validators as val
    if astype is None:
        datatype = megara_inferr_datatype(obj)
        _logger.debug(f'check object as it says it is ({datatype})')
        thistype = datatype
    else:
        _logger.debug(f'check object as {astype}')
        thistype = astype
    checker = val.check_as_datatype(thistype)
    res = checker(obj, level=level)
    return res
