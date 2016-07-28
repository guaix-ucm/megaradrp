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

"""Data model for MEGARA"""

from __future__ import division

import logging
from numina.flow.datamodel import DataModel
from numina.core import DataFrame, ObservationResult


_logger = logging.getLogger('numina.recipes.megara')


class MegaraDataModel(DataModel):
    """Data model of MEGARA.

    Empty for the moment"""

    def __init__(self):
        # Keys
        self._meta = {
            'texp': ('EXPTIME', None),
            'vph': ('VPH', 'undefined'),
            'obsmode': ('OBSMODE', 'undefined'),
            'tstamp': ('TSTAMP', 'undefined')
        }

    def gather_info_dframe(self, img):
        with img.open() as hdulist:
            info = self.gather_info_hdu(hdulist)
        return info

    def gather_info_hdu(self, hdulist):
        values = {}
        values['n_ext'] = len(hdulist)
        extnames = [hdu.header.get('extname', '') for hdu in hdulist[1:]]
        values['name_ext'] = ['PRIMARY'] + extnames
        for key, val in self._meta.items():
            values[key] = hdulist[0].header.get(val[0], val[1])

        return values

    def gather_info(self, recipeinput):
        klass = recipeinput.__class__
        metadata = {}
        for key in klass.stored():
            val = getattr(recipeinput, key)
            if isinstance(val, DataFrame):
                metadata[key] = self.gather_info_dframe(val)
            elif isinstance(val, ObservationResult):
                metas = []
                for f in val.images:
                    metas.append(self.gather_info_dframe(f))
                metadata[key] = metas
            else:
                pass
        return metadata
