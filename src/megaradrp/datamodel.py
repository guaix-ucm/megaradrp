#
# Copyright 2016-2025 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

"""Data model for MEGARA"""

import pkgutil
import logging

import astropy.io.fits as fits
from io import StringIO
from numina.datamodel import DataModel
from numina.keydef import QueryAttribute, KeyDefinition
from numina.util.convert import convert_date

from megaradrp.datatype import MegaraDataType, DataOrigin
import megaradrp.instrument.constants as cons

_logger = logging.getLogger(__name__)


class MegaraDataModel(DataModel):
    """Data model of MEGARA images"""

    query_attrs = {
        "vph": QueryAttribute("vph", str),
        "insmode": QueryAttribute("insmode", str),
        "insconf": QueryAttribute("insconf", str),
        "speclamp": QueryAttribute("speclamp", str),
        "temp": QueryAttribute("temp", float),
        "confid": QueryAttribute("confid", str),
    }

    meta_info_headers = [
        "instrument",
        "object",
        "observation_date",
        "uuid",
        "type",
        "mode",
        "exptime",
        "darktime",
        "insconf",
        "blckuuid",
        "quality_control",
        "vph",
        "insmode",
    ]

    db_info_keys = [
        "instrument",
        "object",
        "observation_date",
        "uuid",
        "type",
        "mode",
        "exptime",
        "darktime",
        "insconf",
        "blckuuid",
        "quality_control",
        "vph",
        "insmode",
    ]

    db_info_keys_extra = ["vph", "insmode"]

    meta_dinfo_headers = [
        "exptime",
        "observation_date",
        "vph",
        "vphpos",
        "insmode",
        "focus",
        "osfilter",
        "uuid",
        "temp",
        "block_uuid",
        "insconf_uuid",
        "speclamp",
        "imgid",
    ]

    PLATESCALE = cons.GTC_FC_A_PLATESCALE.value

    def __init__(self):

        instrument_mappings = {
            "date_obs": ("DATE-OBS", 0, convert_date),
            "insconf": "insconf",
            "blckuuid": "blckuuid",
            "vph": ("VPH", "undefined"),
            "vphpos": ("VPHWHPOS", "undefined"),
            "focus": ("FOCUS", "undefined"),
            "osfilter": ("OSFILTER", "undefined"),
            "temp": ("SENTEMP4", 0.0),
            "speclamp": ("SPECLAMP", "undefined"),
            "confid": KeyDefinition("CONFID", ext="FIBERS"),
        }
        super(MegaraDataModel, self).__init__("MEGARA", instrument_mappings)

    def get_fiberconf(self, img):
        """Obtain FiberConf from image"""
        return get_fiberconf(img)

    def fiber_scale_unit(self, img, unit=False):
        return fiber_scale_unit(img, unit=unit)


def fiber_scale_unit(img, unit=False):
    funit = img["FIBERS"].header.get("FUNIT", "arcsec")

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

    main_insmode = img[0].header.get("INSMODE", "")

    if "FIBERS" in img:
        # We have a 'fibers' extension
        # Information os there
        hdr_fiber = img["FIBERS"].header
        return read_fibers_extension(hdr_fiber, insmode=main_insmode)
    else:
        return get_fiberconf_default(main_insmode)


def create_default_fiber_header(insmode):
    """Obtain default FIBER header"""
    if insmode == "LCB":
        slit_file = "lcb_default_header.txt"
    elif insmode == "MOS":
        slit_file = "mos_default_header.txt"
    else:
        # Read fiber info from headers
        raise ValueError(f"Invalid INSMODE {insmode}")

    data = pkgutil.get_data("megaradrp.instrument.configs", slit_file)
    default_hdr = StringIO(data.decode("utf8"))
    hdr_fiber = fits.header.Header.fromfile(default_hdr)
    return hdr_fiber


def get_fiberconf_default(insmode):
    """Obtain default FiberConf object"""
    hdr_fiber = create_default_fiber_header(insmode)
    return read_fibers_extension(hdr_fiber)


def read_fibers_extension(hdr, insmode="LCB"):
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
    """Description of a MEGARA image: datatype, origin, dates and processing"""
    prim = hdulist[0].header
    instrument = prim.get("INSTRUME", "unknown")
    # date_obs = convert_date(prim.get("DATE-OBS"))
    date_obs = prim.get("DATE-OBS")
    img_uuid = prim.get("UUID")
    insconf = prim.get("INSCONF", "undefined")

    datatype = megara_inferr_datetype_from_image(hdulist)

    obs = {}
    proc = {}
    if datatype.value < MegaraDataType.IMAGE_PROCESSED.value:
        origin = DataOrigin.OBSERVED
    else:
        origin = DataOrigin.PROCESSED

    return {
        "instrument": instrument,
        "datatype": datatype,
        "origin": origin,
        "uuid": img_uuid,
        "insconf": insconf,
        "observation": obs,
        "observation_date": date_obs,
        "processing": proc,
    }


def megara_inferr_datatype(obj):
    """The MegaraDataType of an image or a serialized structured calibration.

    Parameters
    ----------
    obj : astropy.io.fits.HDUList or dict
        An image, see :func:`megara_inferr_datetype_from_image`, or a
        structured calibration serialized as a dictionary, see
        :func:`megara_inferr_datetype_from_dict`.
    """

    if isinstance(obj, fits.HDUList):
        return megara_inferr_datetype_from_image(obj)
    elif isinstance(obj, dict):
        return megara_inferr_datetype_from_dict(obj)
    else:
        raise TypeError(f"I don't know how to inferr datatype from {obj}")


def megara_inferr_datetype_from_dict(obj):
    """The MegaraDataType of a structured calibration, from its type_fqn"""
    # this comes from JSON
    dtype = obj["type_fqn"]
    if dtype in ["megaradrp.products.tracemap.TraceMap"]:
        return MegaraDataType.TRACE_MAP
    elif dtype in ["megaradrp.products.modelmap.ModelMap"]:
        return MegaraDataType.MODEL_MAP
    elif dtype in ["megaradrp.products.wavecalibration.WavelengthCalibration"]:
        return MegaraDataType.WAVE_CALIB
    else:
        return MegaraDataType.UNKNOWN


# The datatypes deduced from the shape of the image
_SHAPE_DATATYPES = {
    (4212, 4196): MegaraDataType.IMAGE_RAW,
    (4112, 4096): MegaraDataType.IMAGE_PROCESSED,
    (623, 4096): MegaraDataType.RSS_PROCESSED,  # LCB
    (644, 4096): MegaraDataType.RSS_PROCESSED,  # MOS
    (623, 4300): MegaraDataType.RSS_WL_PROCESSED,  # LCB
    (644, 4300): MegaraDataType.RSS_WL_PROCESSED,  # MOS
    (4300,): MegaraDataType.SPEC_PROCESSED,
}

# The datatypes of the raw images, from OBSMODE
_OBSMODE_DATATYPES = {
    "MegaraSuccess": MegaraDataType.IMAGE_TEST,
    "MegaraFail": MegaraDataType.IMAGE_TEST,
    "MegaraBiasImage": MegaraDataType.IMAGE_BIAS,
    "MegaraDarkImage": MegaraDataType.IMAGE_DARK,
    "MegaraSlitFlat": MegaraDataType.IMAGE_SLITFLAT,
    "MegaraFiberFlatImage": MegaraDataType.IMAGE_FLAT,
    "MegaraTraceMap": MegaraDataType.IMAGE_FLAT,
    "MegaraModelMap": MegaraDataType.IMAGE_FLAT,
    "MegaraArcCalibration": MegaraDataType.IMAGE_COMP,
    "MegaraTwilightFlatImage": MegaraDataType.IMAGE_TWILIGHT,
    "MegaraLcbImage": MegaraDataType.IMAGE_TARGET,
    "MegaraMosImage": MegaraDataType.IMAGE_TARGET,
    "MegaraMosStdStar": MegaraDataType.IMAGE_TARGET,
    "MegaraExtinctionStar": MegaraDataType.IMAGE_TARGET,
    "MegaraLcbStdStar": MegaraDataType.IMAGE_TARGET,
    "MegaraSensitivityStar": MegaraDataType.IMAGE_TARGET,
    "MegaraFocusTelescope": MegaraDataType.IMAGE_TARGET,
    "MegaraLcbAcquisition": MegaraDataType.IMAGE_TARGET,
    "MegaraMosAcquisition": MegaraDataType.IMAGE_TARGET,
    "MegaraBadPixelMask": MegaraDataType.IMAGE_RAW,
    "MegaraFocusSpectrograph": MegaraDataType.IMAGE_RAW,
}


def _datatype_group(datatype):
    """The datatype that describes the shape of the images of `datatype`"""
    if datatype in [MegaraDataType.MASTER_FLAT, MegaraDataType.MASTER_TWILIGHT]:
        return MegaraDataType.RSS_WL_PROCESSED
    value = datatype.value
    if 100 <= value < 200:
        return MegaraDataType.IMAGE_RAW
    if 200 <= value < 300:
        return MegaraDataType.IMAGE_PROCESSED
    if 300 <= value < 400:
        return MegaraDataType.RSS_PROCESSED
    if 400 <= value < 500:
        return MegaraDataType.RSS_WL_PROCESSED
    if 500 <= value < 600:
        return MegaraDataType.SPEC_PROCESSED
    return MegaraDataType.UNKNOWN


def megara_inferr_datetype_from_image(hdulist):
    """Datatype of a MEGARA image

    The datatype is the one of the keyword IMAGETYP or NUMTYPE, if it
    is a name of MegaraDataType and agrees with the shape of the image.
    Otherwise, the datatype is deduced from the shape of the image and,
    for raw images, from OBSMODE.
    """
    prim = hdulist[0].header
    shape_datatype = _SHAPE_DATATYPES.get(hdulist[0].shape, MegaraDataType.UNKNOWN)

    for key in ["IMAGETYP", "NUMTYPE"]:
        name = prim.get(key)
        if name is None:
            continue
        declared = MegaraDataType.__members__.get(name)
        if declared is None:
            _logger.debug(f"{key}={name} is not a datatype")
        elif shape_datatype == MegaraDataType.UNKNOWN or _datatype_group(declared) == _datatype_group(shape_datatype):
            return declared
        else:
            # The products of a recipe could have the NUMTYPE of its main product
            _logger.debug(f"{key}={name} does not agree with the shape of the image")

    if shape_datatype == MegaraDataType.IMAGE_RAW:
        obsmode = prim.get("OBSMODE", "unknown")
        return _OBSMODE_DATATYPES.get(obsmode, MegaraDataType.UNKNOWN)
    return shape_datatype


def check_obj_megara(obj, astype=None, level=None):
    """Check an image or a structured calibration of MEGARA.

    The object is checked with the checker of its datatype in
    :data:`megaradrp.validators.check_as_datatype`.

    Parameters
    ----------
    obj : astropy.io.fits.HDUList or dict
        The object to check.
    astype : MegaraDataType, optional
        The datatype of the object. If None, it is inferred with
        :func:`megara_inferr_datatype`.
    level : optional
        Not used.

    Raises
    ------
    numina.exceptions.ValidationError
        If the object is not valid.
    """
    import megaradrp.validators as val

    if astype is None:
        datatype = megara_inferr_datatype(obj)
        _logger.debug(f"check object as it says it is ({datatype})")
        thistype = datatype
    else:
        _logger.debug(f"check object as {astype}")
        thistype = astype
    checker = val.check_as_datatype(thistype)
    res = checker(obj, level=level)
    return res
