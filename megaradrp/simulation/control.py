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

import logging

from megaradrp.simulation.actions import megara_sequences
from megaradrp.simulation.factory import PersistentRunCounter

_logger = logging.getLogger("simulation")


class ControlSystem(object):
    """Top level"""
    def __init__(self, factory):
        self._elements = {}

        self.imagecount = PersistentRunCounter('r00%04d.fits')
        self.mode = 'null'
        self.ins = 'MEGARA'
        self.seqs = megara_sequences()
        self.factory = factory
        self.ob_data = dict(count=0, repeat=0, name=None, obsid=1)

    def register(self, name, element):
        self._elements[name] = element

    def get(self, name):
        return self._elements[name]

    def set_mode(self, mode):
        self.mode = mode

    def run(self, exposure, repeat=1):

        if repeat < 1:
            return

        _logger.info('mode is %s', self.mode)
        try:
            thiss = self.seqs[self.mode]
        except KeyError:
            _logger.error('No sequence for mode %s', self.mode)
            raise

        iterf = thiss.run(self, exposure, repeat)

        self.ob_data['repeat'] = repeat
        self.ob_data['name'] = None
        for count, final in enumerate(iterf, 1):
            _logger.info('image %d of %d', count, repeat)
            self.ob_data['name'] = self.imagecount.runstring()
            self.ob_data['count'] = count
            fitsfile = self.factory.create(final, self.ob_data['name'], self)
            _logger.info('save image %s', self.ob_data['name'])
            fitsfile.writeto(self.ob_data['name'], clobber=True)

    def config_info(self):
        return {'ob_data': self.ob_data}

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.imagecount.__exit__(exc_type, exc_val, exc_tb)
