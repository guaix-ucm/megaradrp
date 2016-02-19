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

"""Simple monocromatic simulation"""


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


def simulate_bias_fits(factory, instrument, repeat=1):
    """Simulate a BIAS FITS."""

    # Use instrument and detector!
    det = getattr(instrument, 'detector', None)
    if det:
        detector = det
    else:
        detector = instrument

    for i in range(repeat):
        detector.expose()
        final = detector.readout()
        if det:
            meta = instrument.meta()
        else:
            meta = {'detector': detector.meta()}
        fitsfile = factory.create('bias', meta=meta, data=final)

        yield fitsfile


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
        if det:
            meta = instrument.meta()
        else:
            meta = {'detector': detector.meta()}
        fitsfile = factory.create('dark', meta=meta, data=final)

        yield fitsfile


def simulate_fiber_flat_fits(factory, instrument, wltable, photons_in, exposure, repeat=1):
    """Simulate a FIBER-FLAT"""

    if repeat < 1:
        return

    out = instrument.simulate_focal_plane(wltable, photons_in)

    for i in range(repeat):
        instrument.detector.expose(source=out, time=exposure)

        final = instrument.detector.readout()

        fitsfile = factory.create('fiber-flat', meta=instrument.meta(), data=final)
        yield fitsfile


def simulate_arc_fits(factory, instrument, wltable, photons_in, exposure, repeat=1):
    """Simulate a FLAT"""

    if repeat < 1:
        return

    out = instrument.simulate_focal_plane(wltable, photons_in)

    for i in range(repeat):
        instrument.detector.expose(source=out, time=exposure)

        final = instrument.detector.readout()

        fitsfile = factory.create('arc', meta=instrument.meta(), data=final)
        yield fitsfile


def simulate_focus_fits(factory, instrument, wltable, photons_in, focii, exposure, repeat=1):
    """Simulate a Focus sequence"""

    if repeat < 1:
        return

    for focus in focii:
        instrument.set_focus(focus)

        out = instrument.simulate_focal_plane(wltable, photons_in)

        for i in range(repeat):
            instrument.detector.expose(source=out, time=exposure)

            final = instrument.detector.readout()

            fitsfile = factory.create('focus', meta=instrument.meta(), data=final)
            yield fitsfile
