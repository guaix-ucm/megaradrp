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

def simulate_bias(detector):
    """Simulate a BIAS array."""
    detector.expose(source=0.0, time=0.0)
    final = detector.readout()
    return final

def simulate_dark(detector, exposure):
    """Simulate a DARK array,"""
    detector.expose(source=0.0, time=exposure)
    final = detector.readout()
    return final

def simulate_flat(detector, exposure, source):
    """Simulate a FLAT array,"""
    detector.expose(source=source, time=exposure)
    final = detector.readout()
    return final


def simulate_dark_fits(factory, instrument, exposure, repeat=1):
    """Simulate a DARK FITS."""
    # Use instrument and detector!
    det = getattr(instrument, 'detector', None)
    if det:
        detector = det
    else:
        detector = instrument

    for i in range(repeat):
        detector.expose(time=exposure)
        final = detector.readout()
        fitsfile = factory.create_from_instrument(mode='dark', instrument=detector, data=final, name='dark_%s.fits'%i)

        yield fitsfile


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
        # This is an empty generator
        return iter(())



class MegaraBiasSequence(Sequence):
    def __init__(self):
        super(MegaraBiasSequence, self).__init__('MEGARA', 'bias_image')

    def run(self, control, exposure, repeat):
        instrument = control.get(self.instrument)
        for i in range(repeat):
            instrument.detector.expose()
            final = instrument.detector.readout()
            yield final


class MegaraDarkSequence(Sequence):
    def __init__(self):
        super(MegaraDarkSequence, self).__init__('MEGARA', 'dark_image')

    def run(self, control, exposure, repeat):
        instrument = control.get(self.instrument)

        for i in range(repeat):
            instrument.detector.expose(source=0.0, time=exposure)
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
        if lamp == 'EMPTY':
            raise ValueError('LAMP is %s, exiting' %  lamp)

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
        super(MegaraFiberFlatSequence, self).__init__('fiber_flat_image')

    def lamp_check(self, lamp):
        # TODO: check that this is a cont lamp
        return True


class MegaraTraceMapSequence(MegaraLampSequence):
    def __init__(self):
        super(MegaraTraceMapSequence, self).__init__('trace_map')

    def lamp_check(self, lamp):
        # TODO: check that this is a cont lamp
        return True


class MegaraArcSequence(MegaraLampSequence):
    def __init__(self):
        super(MegaraArcSequence, self).__init__('arc_calibration')

    def lamp_check(self, lamp):
        # TODO: check that this is an arc lamp
        return True


class MegaraSlitFlatSequence(Sequence):
    def __init__(self):
        super(MegaraSlitFlatSequence, self).__init__('MEGARA', 'slit_flat')

    def run(self, control, exposure, repeat):
        # This is an empty generator
        return iter(())


class MegaraTwilightFlatSequence(Sequence):
    def __init__(self):
        super(MegaraTwilightFlatSequence, self).__init__('MEGARA', 'twilight_flat_image')

    def run(self, control, exposure, repeat):
        instrument = control.get(self.instrument)
        telescope = control.get('GTC')
        atm = telescope.inc # Atmosphere model
        # Simulated tw spectrum
        wl_in = instrument.vph.wltable_interp()
        tw_spectrum = atm.twilight_spectrum(wl_in)
        tw_illum = instrument.illumination_in_focal_plane(tw_spectrum)

        out = instrument.simulate_focal_plane(wl_in, tw_illum)
        for i in range(repeat):
            instrument.detector.expose(source=out, time=exposure)
            final = instrument.detector.readout()
            yield final


class MegaraFocusSequence(Sequence):
    def __init__(self):
        super(MegaraFocusSequence, self).__init__('MEGARA', mode='focus_spectrograph')

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
        focii = range(119, 128)

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
    # Keys must correspod to MEGARA ObsMode.key
    seqs['null'] = MegaraNullSequence()
    seqs['bias_image'] = MegaraBiasSequence()
    seqs['dark_image'] = MegaraDarkSequence()
    seqs['fiber_flat_image'] = MegaraFiberFlatSequence()
    seqs['slit_flat'] = MegaraSlitFlatSequence()
    seqs['trace_map'] = MegaraTraceMapSequence()
    seqs['arc_calibration'] = MegaraArcSequence()
    seqs['twilight_flat_image'] = MegaraTwilightFlatSequence()
    seqs['focus_spectrograph'] = MegaraFocusSequence()
    return seqs
