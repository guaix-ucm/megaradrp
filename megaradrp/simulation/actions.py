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
        instrument = control.get(self.instrument)
        cu = control.get('megcalib')

        # Get active lamp
        lamp = cu.current()
        if lamp == 'EMPTY':
            raise ValueError('LAMP is %s, exiting' % lamp)

        # Internal focus
        instrument.set_focus(126.5)
        # FIXME: seting internal value directly
        instrument._internal_focus_factor = 8
        # Simulated arc spectrum

        wl_in, lamp_illum = self.lamp_in_focal_plane(lamp, instrument)
        out = instrument.apply_transmissions(wl_in, lamp_illum)

        for i in range(repeat):
            instrument.detector.expose(source=out, time=exposure)
            final = instrument.detector.readout()
            yield final

    def lamp_in_focal_plane(self, lamp, instrument):
        import numpy as np
        import math

        # print('enter targets in focal plane')
        wl = instrument.vph.wltable_interp()
        # Total number of fibers, defines RSS size
        nfibers = instrument.fibers.nfibers
        # Radius of the fiber (on the focal plane, arcsecs)
        fibrad = instrument.fibers.size
        fibarea = 3 * math.sqrt(3) * fibrad ** 2 / 2.0
        # Result, RSS like (with higher resolution in WL)
        final = np.zeros((nfibers, wl.shape[0]))

        # Centers of profiles
        tab = instrument.focal_plane.get_visible_fibers(instrument.fibers)
        fibid = tab['fibid']
        pos_x = tab['x']
        pos_y = tab['y']
        cover_frac = tab['cover']

        # Extended objects
        # fraction of flux in each fiber
        scales = np.zeros((nfibers,))

        # scales include the cover
        scales[fibid - 1] = cover_frac

        if lamp.illumination:
            scales[fibid - 1] *= lamp.illumination(pos_x, pos_y)

        # FIXME: spectrum is in photons?
        # FIXME: multiply by fiber area
        xx = lamp.flux(wl)
        final += scales[:, np.newaxis] * xx[np.newaxis, :] * fibarea

        return wl, final


class MegaraTwilightFlatSequence(Sequence):
    def __init__(self):
        super(MegaraTwilightFlatSequence, self).__init__('MEGARA', 'twilight_flat_image')

    def run(self, control, exposure, repeat):
        instrument = control.get(self.instrument)
        telescope = control.get('GTC')
        atm = telescope.inc # Atmosphere model
        # Simulated tw spectrum

        wl_in, ns_illum = self.targets_in_focal_plane(control.targets, atm, telescope, instrument)
        out = instrument.apply_transmissions(wl_in, ns_illum)
        for i in range(repeat):
            instrument.detector.expose(source=out, time=exposure)
            final = instrument.detector.readout()
            yield final

    def targets_in_focal_plane(self, targets, atmosphere, telescope, instrument):
        import numpy as np
        import math

        # print('enter targets in focal plane')
        wl = instrument.vph.wltable_interp()
        # Total number of fibers, defines RSS size
        nfibers = instrument.fibers.nfibers
        # Radius of the fiber (on the focal plane, arcsecs)
        fibrad = instrument.fibers.size
        fibarea = 3 * math.sqrt(3) * fibrad ** 2 / 2.0
        # Result, RSS like (with higher resolution in WL)
        final = np.zeros((nfibers, wl.shape[0]))

        # Centers of profiles
        tab = instrument.focal_plane.get_visible_fibers(instrument.fibers)
        fibid = tab['fibid']
        cover_frac = tab['cover']

        # Extended objects
        for _ in [1]:
            # fraction of flux in each fiber
            scales = np.zeros((nfibers,))

            # scales include the cover
            scales[fibid - 1] = cover_frac

            # FIXME, this is hardcoded
            # FIXME: spectrum is in photons?
            # FIXME: multiply by fiber area
            xx = atmosphere.twilight_spectrum(wl)
            final += scales[:, np.newaxis] * xx[np.newaxis, :] * fibarea

        # Multiply by area of the telescope
        final *= telescope.area
        # Multiply by transmission of the optics of the telescope
        final *= telescope.transmission(wl)
        # print('end targets in focal plane')
        return wl, final


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


class MegaraLCBImageSequence(Sequence):
    def __init__(self):
        super(MegaraLCBImageSequence, self).__init__('MEGARA', 'lcb_image')

    def run(self, control, exposure, repeat):
        instrument = control.get(self.instrument)
        telescope = control.get('GTC')
        atm = telescope.inc # Atmosphere model, incoming light
        # Simulated ns spectrum
        # Get targets

        wl_in, ns_illum = self.targets_in_focal_plane(control.targets, atm, telescope, instrument)

        out = instrument.apply_transmissions(wl_in, ns_illum)
        for i in range(repeat):
            instrument.detector.expose(source=out, time=exposure)
            final = instrument.detector.readout()
            yield final

    def targets_in_focal_plane(self, targets, atmosphere, telescope, instrument):
        # simulate_seeing_profile
        from .convolution import rect_c, hex_c, setup_grid
        from scipy.interpolate import RectBivariateSpline
        from scipy import signal
        import numpy as np
        import math

        # print('enter targets in focal plane')
        wl = instrument.vph.wltable_interp()
        # Total number of fibers, defines RSS size
        nfibers = instrument.fibers.nfibers
        # Radius of the fiber (on the focal plane, arcsecs)
        fibrad = instrument.fibers.size
        fibarea = 3 * math.sqrt(3) * fibrad**2 / 2.0
        # Result, RSS like (with higher resolution in WL)
        final = np.zeros((nfibers, wl.shape[0]))

        # Simulation of the fraction of flux in each spaxel
        # By convolution of the seeing profile with the
        # spaxel shape
        # Simulation size
        # FIXME: hardcoded
        # FIXME: this is needed only if there are discrete objects
        xsize = 12.5  # FOR LCB
        ysize = 12.5
        # Pixel size for simulation
        Dx = 0.005
        Dy = 0.005

        xx, yy, xs, ys, xl, yl = setup_grid(xsize, ysize, Dx, Dy)
        # print('generate hex kernel')
        hex_kernel = hex_c(xx, yy, rad=fibrad)

        # print('generate seeing profile')
        seeing_model = atmosphere.seeing
        sc = seeing_model(xx, yy)

        # print('convolve hex kernel with seeing profile')
        # Convolve model gaussian with hex kernel
        convolved = signal.fftconvolve(hex_kernel, sc, mode='same')

        # print('setup 2D interpolator for point-like object')
        # Interpolate
        rbs = RectBivariateSpline(xs, ys, convolved)

        # Centers of profiles
        tab = instrument.focal_plane.get_visible_fibers(instrument.fibers)
        fibid = tab['fibid']
        pos_x = tab['x']
        pos_y = tab['y']
        cover_frac = tab['cover']

        for target in targets:
            # print('object is', target.name)
            center = target.relposition

            # fraction of flux in each fiber
            scales = np.zeros((nfibers,))

            # Offset fiber positions
            offpos0 = pos_x - center[0]
            offpos1 = pos_y - center[1]

            # Check those that are inside the convolved image
            validfibs = rect_c(offpos0, offpos1, 2 * xl, 2 * yl)

            # Sample convolved image
            # In the positions of active fibers
            values = rbs.ev(offpos0[validfibs], offpos1[validfibs]) * Dx * Dy
            # scales include the cover
            scales[fibid[validfibs] - 1] = values * cover_frac[validfibs]

            # If spectrum is in erg s^-1 cm^-2 AA^-1
            # Handle units
            energ_spec = target.spectrum['sed'](wl) * target.profile.factor

            from astropy import units as u
            import astropy.constants as cons

            wl_range = wl * u.AA
            flux = energ_spec * u.erg * u.s ** -1 * u.cm ** -2 * u.AA ** -1

            phot_energy = (cons.h * cons.c / wl_range).to(u.erg)
            nphot_area_t_wl = flux / phot_energy

            nphot_area_t_wl_cgs = nphot_area_t_wl.to(u.cm**-2 * u.s**-1 * u.AA ** -1)
            # FIXME: there is a factor here, we need info from the ETC
            final += scales[:, np.newaxis] * nphot_area_t_wl_cgs[np.newaxis, :] * 1e8

            # print('check, near 1', scales.sum(), final.sum())

        # Extended objects
        for _ in [1]:

            # fraction of flux in each fiber
            scales = np.zeros((nfibers,))

            # scales include the cover
            scales[fibid - 1] = cover_frac

            # FIXME, this is hardcoded
            # FIXME: spectrum is in photons?
            # FIXME: multiply by fiber area
            xx = atmosphere.night_spectrum(wl)
            final += scales[:, np.newaxis] * xx[np.newaxis, :] * fibarea

        # Multiply by area of the telescope
        final *= telescope.area
        # Multiply by transmission of the optics of the telescope
        final *= telescope.transmission(wl)
        # print('end targets in focal plane')
        return wl, final


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
    seqs['lcb_image'] = MegaraLCBImageSequence()
    return seqs
