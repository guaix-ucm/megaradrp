#
# Copyright 2011-2026 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

"""Load image correctors according to present calibrations"""

import warnings

import numina.processing as proc
from numina.processing.getters import get_corrector_for, load_corrector

from megaradrp.processing.trimover import OverscanCorrector, TrimImage
from megaradrp.processing.trimover import GainCorrector
from megaradrp.processing.slitflat import SlitFlatCorrector
from megaradrp.processing.diffuselight import DiffuseLightCorrector

get_corrector_bpm = get_corrector_for("master_bpm", proc.BadPixelCorrector, hwin=0, wwin=3)
get_corrector_bias = get_corrector_for("master_bias", proc.BiasCorrector)
get_corrector_dark = get_corrector_for("master_dark", proc.DarkCorrector)
get_corrector_slit_flat = get_corrector_for("master_slitflat", SlitFlatCorrector)
get_corrector_diffuse_light = get_corrector_for("diffuse_light_image", DiffuseLightCorrector)


def get_corrector_super(rinput, meta, key, correctorclass, datamodel):
    """The corrector of the calibration rinput.<key>

    .. deprecated:: 0.23
        Use numina.processing.getters.load_corrector.
    """
    warnings.warn(
        "get_corrector_super is deprecated, it will be removed in 0.24, "
        "use numina.processing.getters.load_corrector",
        DeprecationWarning,
        stacklevel=2,
    )
    return load_corrector(rinput, key, correctorclass, datamodel)


def get_corrector_overscan(rinput, meta, ins, datamodel):
    detconf = ins.get_property("detector.scan")
    return OverscanCorrector(detconf, datamodel=datamodel, calibid=str(ins.get_device("detector").origin.uuid))


def get_corrector_trimming(rinput, meta, ins, datamodel):
    detconf = ins.get_property("detector.scan")
    return TrimImage(detconf, datamodel=datamodel, calibid=str(ins.get_device("detector").origin.uuid))


def get_corrector_gain(rinput, meta, ins, datamodel):
    """Correct from gain"""
    return GainCorrector(
        detconf=ins.get_property("detector.scan"),
        datamodel=datamodel,
        calibid=str(ins.get_device("detector").origin.uuid),
    )
