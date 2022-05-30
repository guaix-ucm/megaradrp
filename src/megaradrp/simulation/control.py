#
# Copyright 2016-2019 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import logging

from numina.instrument.simulation.factory import PersistentRunCounter

from megaradrp.simulation.actions import megara_sequences
from megaradrp.simulation.efficiency import InterpolFile


_logger = logging.getLogger("simulation")


class PointLike(object):
    def __init__(self, factor=1.0, **kwds):
        self.factor = factor


class Target(object):
    def __init__(self, name, relposition, profile, spectrum):
        self.name = name
        self.relposition = relposition
        self.profile = profile
        self.spectrum = spectrum


def profile_builder(profile_entry):
    profile_type = profile_entry['type']
    profile_args = profile_entry.get('args', {})
    if profile_type == 'point-like':
        return PointLike(**profile_args)
    else:
        raise ValueError("'{0}' not registered", profile_type)


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
        self.targets = None

    def register(self, name, element):
        self._elements[name] = element

    def get(self, name):
        return self._elements[name]

    def set_mode(self, mode):
        self.mode = mode

    def set_targets(self, targets_description):
        # Process targets here

        targets = targets_description['targets']

        tlist = []

        for m in targets:
            entry = targets[m]
            profile_entry = entry['profile']
            profile = profile_builder(profile_entry)
            spectrum_entry = entry['spectrum']
            filename = spectrum_entry['sed']
            factor = spectrum_entry.get('factor', 1.0)
            spectrum_entry['sed'] = InterpolFile(filename, factor=factor)
            spectrum = spectrum_entry
            relposition = entry['relposition']
            tar = Target(m, relposition, profile=profile, spectrum=spectrum)
            tlist.append(tar)

        self.targets = tlist

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
            fitsfile.writeto(self.ob_data['name'], overwrite=True)

    def config_info(self):
        return {'ob_data': self.ob_data}

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.imagecount.__exit__(exc_type, exc_val, exc_tb)
