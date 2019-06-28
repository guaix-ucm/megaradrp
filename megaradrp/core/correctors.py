#
# Copyright 2011-2018 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Load image correctos according to present calibrations"""


import logging


import numina.util.node as node
import numina.processing as proc

from megaradrp.processing.trimover import OverscanCorrector, TrimImage
from megaradrp.processing.trimover import GainCorrector
from megaradrp.processing.slitflat import SlitFlatCorrector


_logger = logging.getLogger(__name__)


def get_corrector_bpm(rinput, meta, ins, datamodel):
    bpm_info = meta.get('master_bpm')
    if bpm_info is not None:
        with rinput.master_bpm.open() as hdul:
            _logger.info('loading BPM')
            mbpm = hdul[0].data
            calibid = datamodel.get_imgid(hdul)
            _logger.debug('BPM image: %s', calibid)
            bpm_corrector = proc.BadPixelCorrector(
                mbpm,
                datamodel=datamodel,
                calibid=calibid
            )
    else:
        _logger.info('BPM not provided, ignored')
        bpm_corrector = node.IdNode()

    return bpm_corrector


def get_corrector_super(rinput, meta, key, correctorclass, datamodel):
    info = meta.get(key)
    req = getattr(rinput, key)
    if req is not None:
        with req.open() as hdul:
            _logger.info('loading %s', key)
            _logger.debug('%s info: %s', key, info)
            datac = hdul['primary'].data
            calibid = datamodel.get_imgid(hdul)
            corrector = correctorclass(datac, datamodel=datamodel,
                                       calibid=calibid)
    else:
        _logger.info('%s not provided, ignored', key)
        corrector = node.IdNode()

    return corrector


def get_corrector_bias(rinput, meta, ins, datamodel):
    key = 'master_bias'
    correctorclass = proc.BiasCorrector
    return get_corrector_super(rinput, meta, key, correctorclass, datamodel)


def get_corrector_dark(rinput, meta, ins, datamodel):
    key = 'master_dark'
    correctorclass = proc.DarkCorrector
    return get_corrector_super(rinput, meta, key, correctorclass, datamodel)


def get_corrector_slit_flat(rinput, meta, ins, datamodel):
    key = 'master_slitflat'
    info = meta.get(key)
    if info is not None:
        req = getattr(rinput, key)
        with req.open() as hdul:
            _logger.info('loading slit flat')
            _logger.debug('%s image: %s', key, info)
            mbpm = hdul[0].data
            calibid = datamodel.get_imgid(hdul)
            corrector = SlitFlatCorrector(mbpm, datamodel, calibid=calibid)
    else:
        _logger.info('%s not provided, ignored', key)
        corrector = node.IdNode()

    return corrector


def get_corrector_overscan(rinput, meta, ins, datamodel):
    detconf = ins.get_property('detector.scan')
    return OverscanCorrector(
        detconf,
        datamodel=datamodel,
        calibid=str(ins.get_device('detector').origin.uuid)
    )


def get_corrector_trimming(rinput, meta, ins, datamodel):
    detconf = ins.get_property('detector.scan')
    return TrimImage(
        detconf,
        datamodel=datamodel,
        calibid=str(ins.get_device('detector').origin.uuid)
    )


def get_corrector_gain(rinput, meta, ins, datamodel):
    """Correct from gain"""
    return GainCorrector(
        detconf=ins.get_property('detector.scan'),
        datamodel=datamodel,
        calibid=str(ins.get_device('detector').origin.uuid)
    )
