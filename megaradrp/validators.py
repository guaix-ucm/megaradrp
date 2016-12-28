#
# Copyright 2016 Universidad Complutense de Madrid
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