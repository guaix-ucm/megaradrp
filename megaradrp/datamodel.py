#
# Copyright 2016-2020 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Data model for MEGARA"""

from __future__ import division

import re
import pkgutil
import enum

import astropy.io.fits as fits
import astropy.table
from six import StringIO
from numina.datamodel import DataModel, QueryAttribute, KeyDefinition
from numina.util.convert import convert_date

import megaradrp.instrument as megins
import megaradrp.instrument.constants as cons

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
        funit = img['FIBERS'].header.get("FUNIT", "arcsec")

        if funit == "arcsec":
            scale = 1
        else:
            scale = self.PLATESCALE
        if unit:
            return (scale, funit)
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
        raise ValueError('Invalid INSMODE {}'.format(insmode))

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
    hdr:
       FITS header
    insmode: str
        default INSMODE

    Returns
    =======
    FibersConf


    """
    import megaradrp.instrument.focalplane as fp
    conf = fp.FibersConf()
    defaults = {}
    defaults['LCB'] = (9, 623)
    defaults['MOS'] = (92, 644)

    if insmode not in ['LCB', 'MOS']:
        raise ValueError('insmode %s not in [LCB, MOS]' % insmode)

    conf.name = hdr.get('INSMODE', insmode)
    conf.conf_id = hdr.get('CONFID', 1)
    conf.nbundles = hdr.get('NBUNDLES', defaults[insmode][0])
    conf.nfibers = hdr.get('NFIBERS', defaults[insmode][1])
    conf.funit = funit = hdr.get("FUNIT", "arcsec")
    # Read bundles

    bun_ids = []
    fib_ids = []
    bundles = conf.bundles
    fibers = conf.fibers

    # loop over everything, count BUN%03d_P and FIB%03d_B
    pattern1 = re.compile(r"BUN(\d+)_P")
    pattern2 = re.compile(r"FIB(\d+)_B")
    for key in hdr:
        bun_re = pattern1.match(key)
        fib_re = pattern2.match(key)
        if bun_re:
            bun_idx = int(bun_re.group(1))
            bun_ids.append(bun_idx)
        elif fib_re:
            fib_idx = int(fib_re.group(1))
            fib_ids.append(fib_idx)

    for i in bun_ids:
        bb = fp.BundleConf()
        bb.id = i
        bb.target_priority = hdr["BUN%03d_P" % i]
        bb.target_name = hdr["BUN%03d_I" % i]
        bb.target_type = fp.TargetType[hdr["BUN%03d_T" % i]]
        bb.enabled = hdr.get("BUN%03d_E" % i, True)
        bb.x = hdr.get("BUN%03d_X" % i, 0.0)
        bb.y = hdr.get("BUN%03d_Y" % i, 0.0)
        bb.pa = hdr.get("BUN%03d_O" % i, 0.0)
        bb.fibers = {}
        bundles[i] = bb

    for fibid in fib_ids:
        ff = fp.FiberConf()
        ff.fibid = fibid

        # Coordinates
        ff.d = hdr["FIB%03d_D" % fibid]
        ff.r = hdr["FIB%03d_R" % fibid]
        ff.o = 0 #hdr["FIB%03d_O" % fibid]
        # Active
        ff.inactive = not hdr["FIB%03d_A" % fibid]

        # Coordinates XY
        ff.x = hdr["FIB%03d_X" % fibid]
        ff.y = hdr["FIB%03d_Y" % fibid]

        ff.bundle_id = hdr["FIB%03d_B" % fibid]
        ff.name = hdr.get("FIB%03d_N" % fibid, 'unknown')

        ff.w1 = hdr.get("FIB%03dW1" % fibid, None)
        ff.w2 = hdr.get("FIB%03dW2" % fibid, None)

        # Validity
        if ff.inactive:
            ff.valid = False
        else:
            ff.valid = hdr.get("FIB%03d_V" % fibid, True)

        bundles[ff.bundle_id].fibers[ff.fibid] = ff
        fibers[ff.fibid] = ff

    return conf


class MegaraDataType(enum.Enum):
    UNKNOWN = 1
    IMAGE_RAW = 100
    IMAGE_BIAS = 102
    IMAGE_DARK = 103
    IMAGE_SLITFLAT = 104
    IMAGE_FLAT = 105
    IMAGE_COMP = 106
    #
    IMAGE_TWILIGHT = 107
    IMAGE_TEST = 109
    IMAGE_TARGET = 150
    #
    IMAGE_PROCESSED = 200
    MASTER_BPM = 201
    MASTER_BIAS = 202
    MasterBias = 202 # Alias
    MASTER_DARK = 203
    MASTER_SLITFLAT = 204
    DIFFUSE_LIGHT = 211
    #
    RSS_PROCESSED = 300
    MASTER_FLAT = 305
    MasterFiberFlat = 305 # Alias
    MASTER_TWILIGHT = 306
    SPEC_PROCESSED = 400
    MASTER_SENSITIVITY = 403
    STRUCT_PROCESSED = 500
    TRACE_MAP = 501
    MODEL_MAP = 502


class DataOrigin(enum.Enum):
    UNKNOWN = 0
    OBSERVED = 1
    PROCESSED = 2
    GENERATED = 3


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
        datatype = megara_inferr_imagetype(hdulist)
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


def megara_inferr_imagetype(hdulist):
    IMAGE_RAW_SHAPE = (4212, 4196)
    IMAGE_PROC_SHAPE = (4112, 4096)
    RSS_IFU_PROC_SHAPE = (4300, 623)
    RSS_MOS_PROC_SHAPE = (4300, 644)
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


def check_3(obj, level=None):
    print('checker3 for {}, level={}'.format(obj, level))


def check_raw(obj, level=None):
    print('check RAW for {}, level={}'.format(obj, level))


_megara_checkers = {}
_megara_checkers[MegaraDataType.UNKNOWN] = check_3
_megara_checkers[MegaraDataType.IMAGE_RAW] = check_raw
_megara_checkers[MegaraDataType.IMAGE_BIAS] = check_raw
_megara_checkers[MegaraDataType.IMAGE_DARK] = check_raw
_megara_checkers[MegaraDataType.IMAGE_SLITFLAT] = check_raw
_megara_checkers[MegaraDataType.IMAGE_FLAT] = check_raw
_megara_checkers[MegaraDataType.IMAGE_COMP] = check_raw
#
_megara_checkers[MegaraDataType.IMAGE_TWILIGHT] = check_raw
_megara_checkers[MegaraDataType.IMAGE_TEST] = check_raw
_megara_checkers[MegaraDataType.IMAGE_TARGET] = check_raw
#
_megara_checkers[MegaraDataType.IMAGE_PROCESSED] = check_3
_megara_checkers[MegaraDataType.MASTER_BPM] = check_3
_megara_checkers[MegaraDataType.MASTER_BIAS] = check_3
_megara_checkers[MegaraDataType.MASTER_DARK] = check_3
_megara_checkers[MegaraDataType.MASTER_SLITFLAT] = check_3
_megara_checkers[MegaraDataType.DIFFUSE_LIGHT] = check_3
#
_megara_checkers[MegaraDataType.RSS_PROCESSED] = check_3
_megara_checkers[MegaraDataType.MASTER_FLAT] = check_3
_megara_checkers[MegaraDataType.MASTER_TWILIGHT] = check_3
_megara_checkers[MegaraDataType.SPEC_PROCESSED] = check_3
_megara_checkers[MegaraDataType.MASTER_SENSITIVITY] = check_3
_megara_checkers[MegaraDataType.STRUCT_PROCESSED] = check_3
_megara_checkers[MegaraDataType.TRACE_MAP] = check_3
_megara_checkers[MegaraDataType.MODEL_MAP] = check_3


def check_as_datatype(datatype):
    return _megara_checkers[datatype]
