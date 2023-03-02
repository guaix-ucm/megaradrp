#
# Copyright 2017-2020 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Enumerated types for instrument configurations of MEGARA"""

import enum


class VphName(enum.Enum):
    VPH405_LR = 1405
    VPH480_LR = 1480
    VPH570_LR = 1570
    VPH675_LR = 1675
    VPH799_LR = 1799
    VPH890_LR = 1890
    VPH410_MR = 2410
    VPH443_MR = 2443
    VPH481_MR = 2481
    VPH521_MR = 2521
    VPH567_MR = 2567
    VPH617_MR = 2617
    VPH656_MR = 2656
    VPH712_MR = 2712
    VPH777_MR = 2777
    VPH926_MR = 2926
    VPH665_HR = 3665
    VPH863_HR = 3863


class BundleType(enum.Enum):
    """Types of bundles"""
    LCB = 1
    RP = 2
    SKY = 3


class FiberPatternType(enum.Enum):
    """Types of fiber pattern on bundles"""
    FIXED = 1
    RP = 2
    SKY = 3


class SlitType(enum.Enum):
    """Types of slits"""
    LCB = 1
    MOS = 2


class SlitPosition(enum.Enum):
    """Positions of the pseudo slit"""
    LCB = 1
    MOS = 2
    OPEN = 3


class TargetType(enum.Enum):
    """Possible targets in a fiber bundle"""
    SOURCE = 1
    UNKNOWN = 2
    UNASSIGNED = 3
    SKY = 4
    REFERENCE = 5
    # aliases for the other fields
    STAR = 5
    BLANK = 4


# Related with Robotic positioner sequences

FIBERMOS_MOV_MAXR1 = 194086
FIBERMOS_MOV_MAXR2 = 0
FIBERMOS_MOV_MINR1 = 0
FIBERMOS_MOV_MINR2 = -15360


class RpMovAxis(enum.Enum):
    RP_MOV_AXIS1 = 1
    RP_MOV_AXIS2 = 2


class SequenceType(enum.Enum):
    SEQUENCE_POS_DEPOS = 1
    SEQUENCE_DEPOS_ONLY = 2
