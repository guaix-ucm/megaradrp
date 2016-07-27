#
# Copyright 2015-2016 Universidad Complutense de Madrid
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

import math
import logging

import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy import signal
from astropy import units as u
import astropy.constants as cons

from .convolution import rect_c, hex_c, setup_grid


_logger = logging.getLogger("megaradrp.simulation")


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

        wl_in, lamp_illum = self.lamp_in_focal_plane(lamp, instrument)
        out = instrument.apply_transmissions(wl_in, lamp_illum)

        for i in range(repeat):
            instrument.detector.expose(source=out, time=exposure)
            final = instrument.detector.readout()
            yield final

    def lamp_check(self, lamp):
        raise NotImplemented

    def lamp_in_focal_plane(self, lamp, instrument):

        # print('enter targets in focal plane')
        wl = instrument.vph.wltable_interp()
        # Total number of fibers, defines RSS size
        nfibers = instrument.fiberset.nfibers
        # Radius of the fiber (on the focal plane, arcsecs)
        fibrad = instrument.fiberset.size
        fibarea = instrument.fiberset.area
        # Result, RSS like (with higher resolution in WL)
        final = np.zeros((nfibers, wl.shape[0]))

        # Centers of profiles
        tab = instrument.focal_plane.get_visible_fibers(instrument.fiberset.name)
        fibid = tab['fibid']
        pos_x = tab['x']
        pos_y = tab['y']
        cover_frac = tab['cover']

        # fraction of flux in each fiber
        scales = np.zeros((nfibers,))

        # scales include the cover
        scales[fibid - 1] = cover_frac

        if lamp.illumination:
            scales[fibid - 1] *= lamp.illumination(pos_x, pos_y)

        photon_energy = (cons.h * cons.c / wl)

        sed = lamp.flux(wl)
        area = 1.0 * u.cm**2 # FIXME, area of the lamp?
        nphot_area_t_wl_cgs = (sed / photon_energy * area * fibarea).to(u.s**-1 * u.micron**-1)

        final += scales[:, np.newaxis] * nphot_area_t_wl_cgs[np.newaxis, :]

        return wl, final


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


class MegaraSlitFlatSequence(MegaraLampSequence):
    def __init__(self):
        super(MegaraSlitFlatSequence, self).__init__('slit_flat')

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


class MegaraTwilightFlatSequence(Sequence):
    def __init__(self):
        super(MegaraTwilightFlatSequence, self).__init__('MEGARA', 'twilight_flat_image')

    def run(self, control, exposure, repeat):
        instrument = control.get(self.instrument)
        telescope = control.get('GTC')
        atm = telescope.inc # Atmosphere model
        # Simulated tw spectrum

        targets2 = [atm.twilight_spectrum]

        wl_in, ns_illum = all_targets_in_focal_plane([], targets2, [], atm, telescope, instrument)

        out1 = instrument.apply_transmissions_only(wl_in, ns_illum)
        out2 = instrument.project_rss(wl_in, out1)

        for i in range(repeat):
            instrument.detector.expose(source=out2, time=exposure)
            final = instrument.detector.readout()
            yield final


class MegaraFocusSequence(MegaraLampSequence):
    def __init__(self):
        super(MegaraFocusSequence, self).__init__(name='focus_spectrograph')

    def run(self, control, exposure, repeat):
        instrument = control.get(self.instrument)
        cu = control.get('megcalib')

        # Get active lamp
        lamp = cu.current()
        self.lamp_check(lamp)
        # Simulated arc spectrum
        wl_in, lamp_illum = self.lamp_in_focal_plane(lamp, instrument)

        # FIXME, hardcoded
        focii = range(119, 128)

        for focus in focii:
            instrument.set_focus(focus)
            out = instrument.apply_transmissions(wl_in, lamp_illum)

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

        # Get targets
        targets1 = control.targets
        targets2 = [atm.night_spectrum]

        wl_in, ns_illum = all_targets_in_focal_plane(targets1, targets2, [], atm, telescope, instrument)

        out1 = instrument.apply_transmissions_only(wl_in, ns_illum)
        out2 = instrument.project_rss(wl_in, out1)

        for i in range(repeat):
            instrument.detector.expose(source=out2, time=exposure)
            final = instrument.detector.readout()
            yield final


class MegaraMOSAcquisitionSequence(Sequence):
    def __init__(self):
        super(MegaraMOSAcquisitionSequence, self).__init__('MEGARA', 'MEGARA_MOS_ACQUISITION')

    def run(self, control, exposure, repeat):
        instrument = control.get(self.instrument)
        telescope = control.get('GTC')

        atm = telescope.inc # Atmosphere model, incoming light

        # Get targets
        targets1 = control.targets
        targets2 = [atm.night_spectrum]

        wl_in, ns_illum = all_targets_in_focal_plane(targets1, targets2, [], atm, telescope, instrument)

        out1 = instrument.apply_transmissions_only(wl_in, ns_illum)
        out2 = instrument.project_rss(wl_in, out1)

        for i in range(repeat):
            instrument.detector.expose(source=out2, time=exposure)
            final = instrument.detector.readout()
            yield final


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
    seqs['MEGARA_MOS_ACQUISITION'] = MegaraMOSAcquisitionSequence()
    return seqs


def simulate_point_like_profile(seeing_model, fibrad, angle=0.0, xsize=12.5, ysize=12.5):
    # Simulation of the fraction of flux in each spaxel
    # By convolution of the seeing profile with the
    # spaxel shape
    # Simulation size
    # xsize, ysize
    # Pixel size for simulation
    Dx = 0.005
    Dy = 0.005

    xx, yy, xs, ys, xl, yl = setup_grid(xsize, ysize, Dx, Dy)
    # print('generate hex kernel')
    # fibrad without units
    hex_kernel = hex_c(xx, yy, rad=fibrad, ang=angle / 180.0 * math.pi)

    # print('generate seeing profile')
    sc = seeing_model(xx, yy)

    # print('convolve hex kernel with seeing profile')
    # Convolve model gaussian with hex kernel
    convolved = signal.fftconvolve(hex_kernel, sc, mode='same')

    # print('setup 2D interpolator for point-like object')
    # Interpolate
    rbs = RectBivariateSpline(xs, ys, convolved)

    def fraction_of_flux(offpos0, offpos1):

        # Check those that are inside the convolved image
        result = np.zeros_like(offpos0)
        validfibs = rect_c(offpos0, offpos1, 2 * xl, 2 * yl)

        # FIXME:
        #import matplotlib.pyplot as plt
        #plt.imshow(sc, extent=[-xl-0.5*Dx, xl+0.5*Dx, -yl-0.5*Dy,yl+0.5*Dy], origin='lower')
        #plt.scatter(offpos0[validfibs], offpos1[validfibs])
        #plt.show()

        # Sample convolved image
        # In the positions of active fibers
        values = rbs.ev(offpos0[validfibs], offpos1[validfibs]) * Dx * Dy
        # scales include the cover
        result[validfibs] = values
        return result

    return fraction_of_flux


def all_targets_in_focal_plane(t1, t2, t3, atmosphere, telescope, instrument):
    # simulate_seeing_profile

    # print('enter targets in focal plane')
    # WL in microns
    wl = instrument.vph.wltable_interp()
    # Total number of fibers, defines RSS size
    nfibers = instrument.fiberset.nfibers
    # Radius of the fiber (on the focal plane, arcsecs)
    fibrad = instrument.fiberset.size
    fibarea = instrument.fiberset.area

    # Result, RSS like (with higher resolution in WL)
    final = np.zeros((nfibers, wl.shape[0]))

    # Centers of profiles
    if instrument.fiberset.name == 'MOS':
        # get fiber mos
        fibermos = instrument.get_device('MEGARA.MOS')
        ids, pos = fibermos.robots_in_focal_plane()
        pos = np.array(pos)

        import scipy.spatial.kdtree as kdtree
        maxdis = 20.0
        base = kdtree.KDTree(pos)
        for target in t1:
            _logger.debug('find robot nearest to %s', target.relposition)
            result = base.query(target.relposition, distance_upper_bound=maxdis)
            if result[1] >= instrument.fiberset.nfibers:
                _logger.debug('no robot nearest to point within %f', maxdis)
                continue
            val = ids[result[1]]
            _logger.debug('robot is number % d', val)
            robot = instrument.get_device('MEGARA.MOS.RoboticPositioner_%d' % val)
            fibid, allpos =  robot.fibers_in_focal_plane()
            tab = instrument.focal_plane.filter_visible_fibers(fibid, allpos)
            fibid = tab['fibid']
            pos_x = tab['x']
            pos_y = tab['y']
            cover_frac = tab['cover']
            subfinal = final[fibid - 1]
            # FIXME: check angles
            rotang = 90 - robot.pa

            add_target_mos(target, subfinal, wl, fibrad, pos_x, pos_y, rotang, telescope, atmosphere)
    else:
        lcb = instrument.get_device('MEGARA.LCB')
        fibid, allpos = lcb.fibers_in_focal_plane()
        tab = instrument.focal_plane.filter_visible_fibers(fibid, allpos)
        fibid = tab['fibid']
        pos_x = tab['x']
        pos_y = tab['y']
        cover_frac = tab['cover']
        subfinal = final[fibid - 1]
        add_targets_lcb(t1, subfinal, wl, fibrad, pos_x, pos_y, telescope, atmosphere)


    add_sky(t2, subfinal, wl, fibarea, telescope)

    # Uniform lamp objects
    add_lamp(t3, subfinal, wl, fibarea, pos_x, pos_y)

    subfinal *= cover_frac[:, np.newaxis]

    final[fibid-1] = subfinal

    return wl, final


def add_targets_lcb(targets, subfinal, wl, fibrad, pos_x, pos_y, telescope, atmosphere):

    fraction_of_flux = simulate_point_like_profile(atmosphere.seeing, fibrad.value)

    # res = np.zeros_like(subfinal)

    airmass = 1.0
    extinction = np.power(10, -0.4 * airmass * atmosphere.extinction(wl))

    energy_unit = u.erg * u.s**-1 * u.cm**-2 * u.AA **-1
    photon_energy = (cons.h * cons.c / wl)

    for target in targets:
        _logger.debug('object is name is %s', target.name)
        center = target.relposition
        # fraction of flux in each fiber
        #scales = np.zeros((nfibers,))

        # Offset fiber positions
        offpos0 = pos_x - center[0]
        offpos1 = pos_y - center[1]

        f_o_f = fraction_of_flux(offpos0, offpos1)
        # If spectrum is in erg s^-1 cm^-2 AA^-1
        # Handle units
        sed = target.spectrum['sed'](wl) * extinction * energy_unit

        nphotons = sed / photon_energy * telescope.area * telescope.transmission(wl)
        subfinal += f_o_f[:, np.newaxis] * nphotons.to(u.s**-1 * u.micron**-1)

    return subfinal


def add_target_mos(target, subfinal, wl, fibrad, pos_x, pos_y, rotang, telescope, atmosphere):

    fraction_of_flux = simulate_point_like_profile(atmosphere.seeing, fibrad.value, angle=rotang, xsize=5.0, ysize=5.0)

    # res = np.zeros_like(subfinal)
    airmass = 1.0
    extinction = np.power(10, -0.4 * airmass * atmosphere.extinction(wl))

    energy_unit = u.erg * u.s**-1 * u.cm**-2 * u.AA **-1
    photon_energy = (cons.h * cons.c / wl)

    _logger.debug('object is %s', target.name)
    center = target.relposition
    # fraction of flux in each fiber
    #scales = np.zeros((nfibers,))

    # Offset fiber positions
    offpos0 = pos_x - center[0]
    offpos1 = pos_y - center[1]

    f_o_f = fraction_of_flux(offpos0, offpos1)
    _logger.debug('fraction of flux recovered %f', f_o_f.sum())
    # If spectrum is in erg s^-1 cm^-2 AA^-1
    # Handle units
    sed = target.spectrum['sed'](wl) * extinction * energy_unit

    nphotons = sed / photon_energy * telescope.area * telescope.transmission(wl)
    subfinal += f_o_f[:, np.newaxis] * nphotons.to(u.s**-1 * u.micron**-1)

    return subfinal

# Uniform atmosphere objects
def add_sky(targets, subfinal, wl, fibarea, telescope):

    energy_unit = u.erg * u.s**-1 * u.cm**-2 * u.AA**-1 * u.arcsec**-2

    photon_energy = (cons.h * cons.c / wl)

    for target in targets:
        sed = target(wl) * energy_unit
        nphotons = sed / photon_energy * fibarea * telescope.area * telescope.transmission(wl)
        subfinal += nphotons.to(u.s**-1 * u.micron**-1)

    return subfinal


# Uniform lamp objects
def add_lamp(lamps, subfinal, wl, fibarea, pos_x, pos_y):

    photon_energy = (cons.h * cons.c / wl)

    for lamp in lamps:
        scales = np.ones_like(pos_x)
        if lamp.illumination:
            scales *= lamp.illumination(pos_x, pos_y)

        sed = lamp.flux(wl)
        area = 1.0 * u.cm ** 2  # FIXME, lamp surface area?
        nphotons = (sed / photon_energy * area * fibarea).to(u.s ** -1 * u.micron ** -1)
        subfinal += scales[:, np.newaxis] * nphotons[np.newaxis, :]

    return subfinal
