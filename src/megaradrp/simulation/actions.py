#
# Copyright 2015-2021 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Sequences for observing modes of MEGARA"""

import math
import logging

import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy import signal
import scipy.spatial
from astropy import units as u
import astropy.constants as cons
from numina.instrument.simulation.actions import Sequence

from .convolution import rect_c, hex_c, setup_grid


_logger = logging.getLogger(__name__)


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
        fitsfile = factory.create_from_instrument(mode='MegaraDarkImage',
                                                  instrument=detector,
                                                  data=final,
                                                  name=f'dark_{i}.fits'
                                                  )

        yield fitsfile


class MegaraSequence(Sequence):
    def __init__(self, mode):
        super(MegaraSequence, self).__init__('MEGARA', mode)


class MegaraNullSequence(Sequence):
    def __init__(self):
        super(MegaraNullSequence, self).__init__('MEGARA', 'null')

    def run(self, control, exposure, repeat):
        # This is an empty generator
        return iter(())


class MegaraBiasSequence(MegaraSequence):
    def __init__(self):
        super(MegaraBiasSequence, self).__init__('MegaraBiasImage')

    def setup_instrument(self, instrument):
        instrument.shutter = 'STOP'

    def run(self, control, exposure, repeat):
        instrument = control.get(self.instrument)

        self.setup_instrument(instrument)

        for i in range(repeat):
            instrument.detector.expose()
            final = instrument.detector.readout()
            yield final


class MegaraDarkSequence(MegaraSequence):
    def __init__(self):
        super(MegaraDarkSequence, self).__init__('MegaraDarkImage')

    def setup_instrument(self, instrument):
        instrument.shutter = 'STOP'

    def run(self, control, exposure, repeat):
        instrument = control.get(self.instrument)

        self.setup_instrument(instrument)

        for i in range(repeat):
            instrument.detector.expose(source=0.0, time=exposure)
            final = instrument.detector.readout()
            yield final


class MegaraLampSequence(Sequence):
    def __init__(self, mode):
        super(MegaraLampSequence, self).__init__('MEGARA', mode)

    def setup_instrument(self, instrument):
        instrument.shutter = 'OPEN'

    def run(self, control, exposure, repeat):
        instrument = control.get(self.instrument)
        cu = control.get('ICM-MEGARA')

        self.setup_instrument(instrument)

        # Get active lamp
        lamp = cu.current()
        if lamp == 'EMPTY':
            raise ValueError(f'LAMP is {lamp}, exiting')

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
        super(MegaraFiberFlatSequence, self).__init__('MegaraFiberFlatImage')

    def lamp_check(self, lamp):
        # TODO: check that this is a cont lamp
        return True


class MegaraTraceMapSequence(MegaraLampSequence):
    def __init__(self):
        super(MegaraTraceMapSequence, self).__init__('MegaraTraceAmp')

    def lamp_check(self, lamp):
        # TODO: check that this is a cont lamp
        return True


class MegaraArcSequence(MegaraLampSequence):
    def __init__(self):
        super(MegaraArcSequence, self).__init__('MegaraArcCalibration')

    def lamp_check(self, lamp):
        # TODO: check that this is an arc lamp
        return True


class MegaraSlitFlatSequence(MegaraLampSequence):
    def __init__(self):
        super(MegaraSlitFlatSequence, self).__init__('MegaraSlitFlat')

    def setup_instrument(self, instrument):
        instrument.shutter = 'OPEN'

        # Internal focus, we need to defocus a lot
        instrument.set_focus(1000)
        # FIXME: seting internal value directly
        instrument._internal_focus_factor = 8

    def run(self, control, exposure, repeat):

        instrument = control.get(self.instrument)
        cu = control.get('ICM-MEGARA')

        self.setup_instrument(instrument)

        # Get active lamp
        lamp = cu.current()
        if lamp == 'EMPTY':
            raise ValueError(f'LAMP is {lamp}, exiting')

        # Simulated arc spectrum
        wl_in, lamp_illum = self.lamp_in_focal_plane(lamp, instrument)
        out = instrument.apply_transmissions(wl_in, lamp_illum)

        for i in range(repeat):
            instrument.detector.expose(source=out, time=exposure)
            final = instrument.detector.readout()
            yield final


class MegaraTwilightFlatSequence(Sequence):
    def __init__(self):
        super(MegaraTwilightFlatSequence, self).__init__('MEGARA', 'MegaraTwilightFlatImage')

    def setup_instrument(self, instrument):
        instrument.shutter = 'OPEN'

    def run(self, control, exposure, repeat):
        instrument = control.get(self.instrument)
        telescope = control.get('GTC')
        atm = telescope.inc # Atmosphere model
        oe = control.get('OE')

        self.setup_instrument(instrument)
        # Simulated tw spectrum

        targets2 = [atm.twilight_spectrum]

        wl_in, ns_illum = all_targets_in_focal_plane([], targets2, [], atm, telescope, oe, instrument)
        out1 = instrument.apply_transmissions_only(wl_in, ns_illum)
        out2 = instrument.project_rss(wl_in, out1)

        for i in range(repeat):
            instrument.detector.expose(source=out2, time=exposure)
            final = instrument.detector.readout()
            yield final


class MegaraFocusSequence(MegaraLampSequence):
    def __init__(self):
        super(MegaraFocusSequence, self).__init__(mode='MegaraFocusSpectrograph')

    def run(self, control, exposure, repeat):
        instrument = control.get(self.instrument)
        instrument.shutter = 'OPEN'
        cu = control.get('ICM-MEGARA')

        # Get active lamp
        lamp = cu.current()
        self.lamp_check(lamp)
        # Simulated arc spectrum
        wl_in, lamp_illum = self.lamp_in_focal_plane(lamp, instrument)

        # FIXME, hardcoded
        focii = range(-1000, 4000, 1000)
        for focus in focii:
            instrument.set_focus(focus)
            out = instrument.apply_transmissions(wl_in, lamp_illum)

            for i in range(repeat):
                instrument.detector.expose(source=out, time=exposure)
                final = instrument.detector.readout()
                yield final

    def lamp_check(self, lamp):
        return True


class MegaraSkyImageSequence(MegaraSequence):
    def __init__(self, mode):
        super(MegaraSkyImageSequence, self).__init__(mode)

    def setup_instrument(self, instrument):
        # Setup instrument
        # Fixed by mode
        instrument.shutter = 'OPEN'

    def run(self, control, exposure, repeat):

        instrument = control.get(self.instrument)

        self.setup_instrument(instrument)

        telescope = control.get('GTC')
        oeng = control.get('OE')
        atm = telescope.inc # Atmosphere model, incoming light

        # Get targets
        targets1 = control.targets
        targets2 = [atm.night_spectrum]

        wl_in, ns_illum = all_targets_in_focal_plane(targets1, targets2, [], atm, telescope, oeng, instrument)

        out1 = instrument.apply_transmissions_only(wl_in, ns_illum)
        out2 = instrument.project_rss(wl_in, out1)

        for i in range(repeat):
            instrument.detector.expose(source=out2, time=exposure)
            final = instrument.detector.readout()
            yield final


class MegaraSkyLCBImageSequence(MegaraSkyImageSequence):
    def __init__(self, mode):
        super(MegaraSkyLCBImageSequence, self).__init__(mode)

    def setup_instrument(self, instrument):
        super(MegaraSkyLCBImageSequence, self).setup_instrument(instrument)
        instrument.insmode = 'LCB'


class MegaraLCBImageSequence(MegaraSkyLCBImageSequence):
    def __init__(self):
        super(MegaraLCBImageSequence, self).__init__('MegaraLcbImage')


class MegaraLCBAcquisitionSequence(MegaraSkyLCBImageSequence):
    def __init__(self):
        super(MegaraLCBAcquisitionSequence, self).__init__('MegaraLcbAcquisition')


class MegaraFocusTelescopeSequence(MegaraSkyLCBImageSequence):
    def __init__(self):
        super(MegaraFocusTelescopeSequence, self).__init__('MegaraFocusTelescope')

    def run(self, control, exposure, repeat):
        telescope = control.get('GTC')

        # FIXME, hardcoded
        focii = range(-3000, 4000, 1000)
        for focus in focii:
            telescope.set_focus(focus)

            for img in super(MegaraFocusTelescopeSequence, self).run(control, exposure, repeat):
                yield img


class MegaraSkyMOSImageSequence(MegaraSkyImageSequence):
    def __init__(self, mode):
        super(MegaraSkyMOSImageSequence, self).__init__(mode)

    def setup_instrument(self, instrument):
        super(MegaraSkyMOSImageSequence, self).setup_instrument(instrument)
        instrument.insmode = 'MOS'


class MegaraMOSImageSequence(MegaraSkyMOSImageSequence):
    def __init__(self):
        super(MegaraMOSImageSequence, self).__init__('MegaraMosImage')


class MegaraMOSAcquisitionSequence(MegaraSkyMOSImageSequence):
    def __init__(self):
        super(MegaraMOSAcquisitionSequence, self).__init__('MegaraMosAcquisition')


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
    seqs['focus_telescope'] = MegaraFocusTelescopeSequence()
    seqs['lcb_image'] = MegaraLCBImageSequence()
    seqs['mos_image'] = MegaraMOSImageSequence()
    seqs['mos_acquisition'] = MegaraMOSAcquisitionSequence()
    seqs['lcb_acquisition'] = MegaraLCBAcquisitionSequence()
    return seqs


def simulate_point_like_profile(seeing_profile, psf, fibrad, angle=0.0, xsize=12.5, ysize=12.5):
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
    sc = seeing_profile(xx, yy)

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

        # import matplotlib.pyplot as plt
        # plt.imshow(sc, extent=[-xl-0.5*Dx, xl+0.5*Dx, -yl-0.5*Dy,yl+0.5*Dy], origin='lower')
        # plt.scatter(offpos0[validfibs], offpos1[validfibs])
        # plt.show()

        # Sample convolved image
        # In the positions of active fibers
        #
        values = rbs.ev(offpos0[validfibs], offpos1[validfibs]) * Dx * Dy
        # scales include the cover
        result[validfibs] = values
        return result

    return fraction_of_flux


def all_targets_in_focal_plane(t1, t2, t3, atmosphere, telescope, oe, instrument):
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
    # zenith_distance = 30 * u.deg
    # airmass = 1.0 / math.cos(zenith_distance.to(u.rad).value)

    if instrument.fiberset.name == 'MOS':
        # get fiber mos
        fibermos = instrument.get_device('MEGARA.MOS')
        ids, pos = fibermos.robots_in_focal_plane()
        pos = np.array(pos)

        maxdis = 20.0
        base = scipy.spatial.KDTree(pos)
        for target in t1:
            _logger.debug('simulate MOS')
            _logger.debug('find robot nearest to %s', target.relposition)
            result = base.query(target.relposition, distance_upper_bound=maxdis)
            if result[1] >= instrument.fiberset.nfibers:
                _logger.debug('no robot nearest to point within %f', maxdis)
                continue
            val = ids[result[1]]
            _logger.debug('robot is number % d', val)
            robot = instrument.get_device('MEGARA.MOS.RoboticPositioner_%d' % val)
            fibid1, allpos1 =  robot.fibers_in_focal_plane()
            tab1 = instrument.focal_plane.filter_visible_fibers(fibid1, allpos1)
            fibid1 = tab1['fibid']
            pos_x = tab1['x']
            pos_y = tab1['y']
            subfinal = final[fibid1 - 1]
            # FIXME: check angles
            rotang = 90 - robot.pa
            subfinal = add_target_mos(target, subfinal, wl, fibrad, pos_x, pos_y, rotang, oe, telescope, atmosphere)
            final[fibid1 - 1] = subfinal

        # Recompute all visible fibers
        fibid, allpos = fibermos.fibers_in_focal_plane()
        tab = instrument.focal_plane.filter_visible_fibers(fibid, allpos)
        fibid = tab['fibid']
        pos_x = tab['x']
        pos_y = tab['y']
        cover_frac = tab['cover']
        subfinal = final[fibid - 1]
    else:
        _logger.debug('simulate LCB')
        lcb = instrument.get_device('MEGARA.LCB')
        fibid, allpos = lcb.fibers_in_focal_plane()
        tab = instrument.focal_plane.filter_visible_fibers(fibid, allpos)
        fibid = tab['fibid']
        pos_x = tab['x']
        pos_y = tab['y']
        cover_frac = tab['cover']
        subfinal = final[fibid - 1]
        add_targets_lcb2(t1, subfinal, wl, fibrad, pos_x, pos_y, oe, telescope, atmosphere)

    add_sky(t2, subfinal, wl, fibarea, telescope)

    # Uniform lamp objects
    add_lamp(t3, subfinal, wl, fibarea, pos_x, pos_y)

    subfinal *= cover_frac[:, np.newaxis]

    final[fibid-1] = subfinal

    return wl, final


def add_targets_lcb(targets, subfinal, wl, fibrad, pos_x, pos_y, oe, telescope, atmosphere):

    airmass = oe.airmass
    extinction = np.power(10, -0.4 * airmass * atmosphere.extinction(wl))
    seeing_fwhm = atmosphere.seeing.fwhm(wl[0], oe.zenith_distance)
    _logger.debug('seeing FWHM is %', seeing_fwhm)
    # This is not really correct, but is simpler
    # than doing the whole PSF convolution
    _logger.debug('internal focus factor %s', telescope.focus_actuator.internal_focus_factor)
    # Scale with focus
    seeing_fwhm_focus = seeing_fwhm * telescope.focus_actuator.internal_focus_factor
    telescope.focus_actuator.focus = 1200
    print('intnl focus factor', telescope.focus_actuator.internal_focus_factor)
    seeing_profile = atmosphere.seeing.profile(seeing_fwhm_focus)
    psf = None
    fraction_of_flux = simulate_point_like_profile(seeing_profile, psf, fibrad.value)

    energy_unit = u.erg * u.s**-1 * u.cm**-2 * u.AA **-1
    photon_energy = (cons.h * cons.c / wl)

    for target in targets:
        _logger.debug('object name is %s', target.name)
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


def add_targets_lcb2(targets, subfinal, wl, fibrad, pos_x, pos_y, oe, telescope, atmosphere):
    """This correction includes DAR"""
    ref_wl = 0.5 * (wl.max() + wl.min())
    _logger.debug('reference wl is %s', ref_wl)
    seeing_fwhm = atmosphere.seeing.fwhm(ref_wl, oe.zenith_distance)
    _logger.debug('seeing FWHM is %s', seeing_fwhm)
    # Scale with focus
    # telescope.focus_actuator.focus = 3000
    _logger.debug('internal focus factor %s', telescope.focus_actuator.internal_focus_factor)
    seeing_fwhm_focus = seeing_fwhm * telescope.focus_actuator.internal_focus_factor
    #
    seeing_profile = atmosphere.seeing.profile(seeing_fwhm_focus)
    psf = None
    fraction_of_flux = simulate_point_like_profile(seeing_profile, psf, fibrad.value)

    airmass = oe.airmass
    zenith_distance = oe.zenith_distance
    extinction = np.power(10, -0.4 * airmass * atmosphere.extinction(wl))

    dar = atmosphere.refraction(zenith_distance, wl, ref_wl).to(u.arcsec).value
    _logger.debug('extreme DAR %s %s', dar[0], dar[-1])

    energy_unit = u.erg * u.s**-1 * u.cm**-2 * u.AA **-1
    photon_energy = (cons.h * cons.c / wl)

    for target in targets:
        _logger.debug('object is %s', target.name)
        center = target.relposition

        # Offset fiber positions
        # Assume DAR affects only one coordinate, this depends on the several
        # different angles
        center_dar = np.zeros((2, wl.size))
        center_dar[0] = center[0] + 0.0 * dar
        center_dar[1] = center[1] + 1.0 * dar
        offpos1_dar = pos_y[:, None] - center_dar[1]
        offpos0_dar = pos_x[:, None] - center_dar[0]

        offpos0 = pos_x - center[0]
        offpos1 = pos_y - center[1]
        f_o_f_dar = fraction_of_flux(offpos0_dar, offpos1_dar)
        f_o_f_wl = f_o_f_dar.sum(axis=0)
        # print f_o_f_dar
        f_o_f = fraction_of_flux(offpos0, offpos1)
        _logger.debug('fraction of flux recovered (no DAR) %f', f_o_f.sum())
        _logger.debug('fraction of flux recovered min %f', f_o_f_wl.min())
        _logger.debug('fraction of flux recovered max %f', f_o_f_wl.max())

        # If spectrum is in erg s^-1 cm^-2 AA^-1
        # Handle units
        sed = target.spectrum['sed'](wl) * extinction * energy_unit

        nphotons = sed / photon_energy * telescope.area * telescope.transmission(wl)
        flux_dar = f_o_f_dar * nphotons.to(u.s**-1 * u.micron**-1)
        # flux_no_dar = f_o_f[:, np.newaxis] * nphotons.to(u.s ** -1 * u.micron ** -1)
        subfinal += flux_dar

    return subfinal


def add_target_mos(target, subfinal, wl, fibrad, pos_x, pos_y, rotang, oe, telescope, atmosphere):

    ref_wl = 0.5 * (wl.max() + wl.min())
    _logger.debug('reference wl is %s', ref_wl)

    seeing_fwhm = atmosphere.seeing.fwhm(ref_wl, oe.zenith_distance)
    _logger.debug('seeing FWHM is %s', seeing_fwhm)
    # telescope.focus_actuator.focus = 3000
    _logger.debug('internal focus factor %s', telescope.focus_actuator.internal_focus_factor)
    seeing_fwhm_focus = seeing_fwhm * telescope.focus_actuator.internal_focus_factor
    seeing_profile = atmosphere.seeing.profile(seeing_fwhm_focus)
    psf = None
    fraction_of_flux = simulate_point_like_profile(seeing_profile, psf, fibrad.value, angle=rotang, xsize=5.0, ysize=5.0)

    airmass = oe.airmass
    _logger.debug('airmass is %s', airmass)
    zenith_distance = oe.zenith_distance
    _logger.debug('Z dis is %s', zenith_distance)

    extinction = np.power(10, -0.4 * airmass * atmosphere.extinction(wl))

    energy_unit = u.erg * u.s**-1 * u.cm**-2 * u.AA **-1
    photon_energy = (cons.h * cons.c / wl)

    _logger.debug('object is %s', target.name)
    center = target.relposition

    ref_wl = 0.5 * (wl.max()  + wl.min())

    _logger.debug('reference wl is %s', ref_wl)
    dar = atmosphere.refraction(zenith_distance, wl, ref_wl).to(u.arcsec).value
    _logger.debug('extreme DAR %s %s', dar[0], dar[-1])
    # Offset fiber positions
    # Assume DAR affects only one coordinate, this depends on the several
    # different angles
    center_dar = np.zeros((2, wl.size))
    center_dar[0] = center[0] + 0.0 * dar
    center_dar[1] = center[1] + 1.0 * dar
    offpos1_dar = pos_y[:, None] - center_dar[1]
    offpos0_dar = pos_x[:, None] - center_dar[0]

    offpos0 = pos_x - center[0]
    offpos1 = pos_y - center[1]
    f_o_f_dar = fraction_of_flux(offpos0_dar, offpos1_dar)
    f_o_f_wl = f_o_f_dar.sum(axis=0)
    # print f_o_f_dar
    f_o_f = fraction_of_flux(offpos0, offpos1)
    _logger.debug('fraction of flux recovered (no DAR) %f', f_o_f.sum())
    _logger.debug('fraction of flux recovered min %f', f_o_f_wl.min())
    _logger.debug('fraction of flux recovered max %f', f_o_f_wl.max())

    # If spectrum is in erg s^-1 cm^-2 AA^-1
    # Handle units
    sed = target.spectrum['sed'](wl) * extinction * energy_unit

    nphotons = sed / photon_energy * telescope.area * telescope.transmission(wl)
    flux_dar = f_o_f_dar * nphotons.to(u.s**-1 * u.micron**-1)
    # flux_no_dar = f_o_f[:, np.newaxis] * nphotons.to(u.s ** -1 * u.micron ** -1)
    subfinal += flux_dar

    return subfinal


# Uniform atmosphere objects
def add_sky(targets, subfinal, wl, fibarea, telescope):

    energy_unit = u.erg * u.s**-1 * u.cm**-2 * u.AA**-1 * u.arcsec**-2

    photon_energy = (cons.h * cons.c / wl)
    _logger.debug('add sky targets')
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
