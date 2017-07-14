#
# Copyright 2011-2017 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# Megara DRP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Megara DRP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Megara DRP.  If not, see <http://www.gnu.org/licenses/>.
#


"""Load image correctos according to present calibrations"""


import logging

from numina.flow.node import IdNode
from numina.flow.processing import BiasCorrector, BadPixelCorrector
from numina.flow.processing import DarkCorrector

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
            bpm_corrector = BadPixelCorrector(mbpm,
                                              datamodel=datamodel,
                                              calibid=calibid)
    else:
        _logger.info('BPM not provided, ignored')
        bpm_corrector = IdNode()

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
        corrector = IdNode()

    return corrector


def get_corrector_bias(rinput, meta, ins, datamodel):
    key = 'master_bias'
    correctorclass = BiasCorrector
    return get_corrector_super(rinput, meta, key, correctorclass, datamodel)


def get_corrector_dark(rinput, meta, ins, datamodel):
    key = 'master_dark'
    correctorclass = DarkCorrector
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
        corrector = IdNode()

    return corrector


def get_corrector_overscan(rinput, meta, ins, datamodel):
    detconf = ins.get('detector.scan')
    return OverscanCorrector(
        detconf,
        datamodel=datamodel,
        calibid=ins.components['detector'].uuid
    )


def get_corrector_trimming(rinput, meta, ins, datamodel):
    detconf = ins.get('detector.scan')
    return TrimImage(
        detconf,
        datamodel=datamodel,
        calibid=ins.components['detector'].uuid
    )


def get_corrector_gain(rinput, meta, ins, datamodel):
    """Correct from gain"""
    return GainCorrector(
        detconf=ins.get('detector.scan'),
        datamodel=datamodel,
        calibid=ins.components['detector'].uuid
    )

