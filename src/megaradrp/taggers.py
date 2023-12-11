#
# Copyright 2015 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

from numina.core.taggers import get_tags_from_full_ob


def tagger_empty(obsres):
    return {}


def tagger_vph(obsres):
    return get_tags_from_full_ob(obsres, reqtags=['vph'])


def tagger_base_image(obsres):
    return get_tags_from_full_ob(obsres, reqtags=['vph', 'insmode'])
