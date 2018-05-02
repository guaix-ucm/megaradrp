#
# Copyright 2016-2017 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Validators for Observing modes"""

import sys

import six

from numina.exceptions import ValidationError


def validate_focus(obresult):
    """Validate FOCUS_SPECTROGRAPH"""
    image_groups = {}
    for idx, frame in enumerate(obresult.frames):
        with frame.open() as img:
            try:
                focus_val = img[0].header['FOCUS']
                # FIXME: This should be done using an scheme for MEGARA
                if not isinstance(focus_val, (int, float)):
                    raise ValidationError("FOCUS must be integer, not {}".format(type(focus_val)))
                if focus_val not in image_groups:
                    image_groups[focus_val] = []
                image_groups[focus_val].append(frame)
            except Exception:
                _type, exc, tb = sys.exc_info()
                six.reraise(ValidationError, exc, tb)

    if len(image_groups) < 2:
        raise ValidationError('We have only {} different focus in OB'.format(len(image_groups)))

    return True