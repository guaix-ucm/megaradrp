#
# Copyright 2015 Universidad Complutense de Madrid
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

"""Sequences for observing modes of MEGARA"""


class Sequence(object):
    def __init__(self, instrument, mode):
        self.instrument = instrument
        self.mode = mode

    def run(self, **kwds):
        raise NotImplemented


class MegaraNullSequence(Sequence):
    def __init__(self):
        super(MegaraNullSequence, self).__init__('MEGARA', 'null')

    def run(self, control, exposure, repeat):
        pass


class MegaraBiasSequence(Sequence):
    def __init__(self):
        super(MegaraBiasSequence, self).__init__('MEGARA', 'bias')

    def run(self, control, exposure, repeat):
        instrument = control.get(self.instrument)
        for i in range(repeat):
            instrument.detector.expose()
            final = instrument.detector.readout()
            yield final


class MegaraDarkSequence(Sequence):
    def __init__(self):
        super(MegaraDarkSequence, self).__init__('MEGARA', 'dark')

    def run(self, control, exposure, repeat):
        instrument = control.get(self.instrument)

        for i in range(repeat):
            instrument.detector.expose(exposure)
            final = instrument.detector.readout()
            yield final


class MegaraLampSequence(Sequence):
    def __init__(self, name):
        super(MegaraLampSequence, self).__init__('MEGARA', name)

    def run(self, control, exposure, repeat):
        instrument = control.get(self.instrument)
        cu = control.get('megcalib')

        # Get active lamp
        lamp = cu.current()
        self.lamp_check(lamp)
        # Simulated arc spectrum
        wl_in = instrument.vph.wltable_interp()
        lamp_illum = instrument.illumination_in_focal_plane(lamp.flux(wl_in), lamp.illumination)

        out = instrument.simulate_focal_plane(wl_in, lamp_illum)
        for i in range(repeat):
            instrument.detector.expose(source=out, time=exposure)
            final = instrument.detector.readout()
            yield final

    def lamp_check(self, lamp):
        raise NotImplemented


class MegaraFiberFlatSequence(MegaraLampSequence):
    def __init__(self):
        super(MegaraFiberFlatSequence, self).__init__('fiberflat')

    def lamp_check(self, lamp):
        # TODO: check that this is a cont lamp
        return True


class MegaraArcSequence(MegaraLampSequence):
    def __init__(self):
        super(MegaraArcSequence, self).__init__('arc')

    def lamp_check(self, lamp):
        # TODO: check that this is an arc lamp
        return True


class MegaraTwilightFlatSequence(Sequence):
    def __init__(self):
        super(MegaraTwilightFlatSequence, self).__init__('MEGARA', 'twilightflat')

    def run(self, control, exposure, repeat):
        instrument = control.get(self.instrument)
        telescope = control.get('GTC')
        atm = telescope.inc # Atmosphere model
        # Simulated tw spectrum
        wl_in = instrument.vph.wltable_interp()
        tw_spectrum = 1e4 * atm.twightlight_spectrum(wl_in)
        tw_illum = instrument.illumination_in_focal_plane(tw_spectrum)

        out = instrument.simulate_focal_plane(wl_in, tw_illum)
        for i in range(repeat):
            instrument.detector.expose(source=out, time=exposure)
            final = instrument.detector.readout()
            yield final


class MegaraFocusSequence(Sequence):
    def __init__(self):
        super(MegaraFocusSequence, self).__init__('MEGARA', mode='focus')

    def run(self, control, exposure, repeat):
        instrument = control.get(self.instrument)
        cu = control.get('megcalib')

        # Get active lamp
        lamp = cu.current()
        self.lamp_check(lamp)
        # Simulated arc spectrum
        wl_in = instrument.vph.wltable_interp()
        lamp_illum = instrument.illumination_in_focal_plane(lamp.flux(wl_in), lamp.illumination)
        # FIXME, hardcoded
        focii = range(100, 200)

        for focus in focii:
            instrument.set_focus(focus)
            out = instrument.simulate_focal_plane(wl_in, lamp_illum)

            for i in range(repeat):
                instrument.detector.expose(source=out, time=exposure)
                final = instrument.detector.readout()
                yield final

    def lamp_check(self, lamp):
        return True


def megara_sequences():
    seqs = {}
    seqs['null'] = MegaraNullSequence()
    seqs['bias'] = MegaraBiasSequence()
    seqs['dark'] = MegaraDarkSequence()
    seqs['fiberflat'] = MegaraFiberFlatSequence()
    seqs['arc'] = MegaraArcSequence()
    seqs['twilightflat'] = MegaraTwilightFlatSequence()
    seqs['focus'] = MegaraFocusSequence()
    return seqs
