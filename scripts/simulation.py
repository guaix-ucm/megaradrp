from __future__ import print_function

import os
import logging
import datetime

import numpy as np
from astropy import units as u
from astropy.io import fits

from megaradrp.simulation.instrument import MegaraInstrument
from megaradrp.simulation.efficiency import EfficiencyFile
from megaradrp.simulation.efficiency import InterpolFile
from megaradrp.simulation.instrument import InternalOptics
from megaradrp.simulation.wheel import VPHWheel
from megaradrp.simulation import lamps
from megaradrp.simulation import calibrationunit
from megaradrp.simulation.actions import megara_sequences
from megaradrp.simulation.telescope import Telescope
from megaradrp.simulation.fiberbundle import FiberBundle
from megaradrp.simulation.lightfiber import LightFiber, FiberSet
from megaradrp.simulation.psslit import PseudoSlit, PseudoSlitSelector
from megaradrp.simulation.focalplane import FocalPlane
from megaradrp.simulation.detector import ReadParams, MegaraDetectorSat
from megaradrp.simulation.vph import MegaraVPH
from megaradrp.simulation.shutter import MegaraShutter
from megaradrp.simulation.fibermos import FiberMOS, RoboticPositioner, LargeCompactBundle
from megaradrp.simulation.cover import MegaraCover

_logger = logging.getLogger("megaradrp.simulation")


# create detector from data
def create_detector():
    _logger.info('create detector')
    DSHAPE = (2056 * 2, 2048 * 2)
    PSCAN = 50
    OSCAN = 50
    #qe = 1.0 * np.ones(DSHAPE)
    qe = fits.getdata('v03/base_qe.fits')
    dcurrent = 3.0 / 3600

    readpars1 = ReadParams(gain=1.0, ron=2.0, bias=1000.0)
    readpars2 = ReadParams(gain=1.0, ron=2.0, bias=1005.0)

    qe_wl = EfficiencyFile('v03/tccdbroad_0.1aa.dat')

    detector = MegaraDetectorSat('Detector',
                                 DSHAPE, OSCAN, PSCAN, qe=qe, qe_wl=qe_wl, dark=dcurrent,
                                 readpars1=readpars1, readpars2=readpars2, bins='11'
                                 )
    return detector


def create_lcb(insconf):

    _logger.info('create LCB')
    layouttable = np.loadtxt('v03/LCB_spaxel_centers.txt')

    fiberset = FiberSet(name='LCB', size=0.31 * u.arcsec, fwhm=3.6)
    # FIXME: a trans object per fiber is very slow
    trans = EfficiencyFile('v03/tfiber_0.1aa_20m.dat')
    for line in layouttable:
        idx =  int(line[3])
        if idx not in fiberset.bundles:
            name = 'FiberBundle_{}'.format(idx)
            fiberset.bundles[idx] = FiberBundle(name, bid=idx)
        fibid = int(line[4])
        name = 'LightFiber_{}'.format(fibid)

        lf = LightFiber(name, fibid, transmission=trans)
        fiberset.bundles[idx].add_light_fiber(lf)
        fiberset.fibers[fibid] = lf

    lcb_pos = {}
    spaxels = insconf.get('LCB.spaxels')
    for entry in spaxels:
        idx = entry[2]
        lcb_pos[idx] = (entry[0], entry[1])
    lcb = LargeCompactBundle('LCB', fiberset, lcb_pos)

    pseudo_slit_lcb = PseudoSlit(name="LCB", insmode='LCB')
    pseudo_slit_lcb.connect_fibers(fiberset, layouttable[:,2])

    return lcb, pseudo_slit_lcb


def create_mos():
    _logger.info('create MOS')
    layouttable = np.loadtxt('v03/MOS_spaxel_centers.txt')
    # Create fiber bundles and light fibers

    fiberset = FiberSet(name='MOS', size=0.31 * u.arcsec, fwhm=3.6)
    # FIXME: a trans object per fiber is very slow
    trans = EfficiencyFile('v03/tfiber_0.1aa_20m.dat')
    for line in layouttable:
        idx =  int(line[3])
        if idx not in fiberset.bundles:
            name = 'FiberBundle_{}'.format(idx)
            fiberset.bundles[idx] = FiberBundle(name, bid=idx, static=False)
        fibid = int(line[4])
        name = 'LightFiber_{}'.format(fibid)

        lf = LightFiber(name, fibid, transmission=trans)
        fiberset.bundles[idx].add_light_fiber(lf)
        fiberset.fibers[fibid] = lf

    # Center of bundles
    fiber_mos = FiberMOS('MOS', fiberset)
    for line in layouttable[3::7]:
        idx =  int(line[3])
        name = 'RoboticPositioner_{}'.format(idx)
        rb = RoboticPositioner(name, id=idx, pos=(line[0], line[1], 0.0), parent=fiber_mos)
        rb.connect_bundle(fiberset.bundles[idx])

    pseudo_slit_mos = PseudoSlit(name="MOS", insmode='MOS')
    pseudo_slit_mos.connect_fibers(fiberset, layouttable[:, 2])

    return fiber_mos, pseudo_slit_mos


def create_wheel():
    _logger.info('create wheel')
    wheel = VPHWheel(capacity=3, name='Wheel')
    _logger.info('create vphs')


    vph_conf = {'wl_range': [3653.0, 4051.0, 4386.0]}
    vph = create_vph_by_data('VPH405_LR',
                             'LR-U',
                             'v03/VPH405_LR2-extra.dat',
                             'v03/VPH405_LR_res.dat',
                             'v03/tvph_0.1aa.dat',
                             conf=vph_conf
                             )
    wheel.put_in_pos(vph, 0)

    vph_conf = {'wl_range': [8800.0, 9262.0, 9686.0]}
    vph = create_vph_by_data('VPH926_MR',
                             'MR-Z',
                             'v03/VPH926_MR.txt',
                             'v03/VPH926_MR_res.dat',
                             'v03/tvph_0.1aa.dat',
                             conf=vph_conf
                             )
    wheel.put_in_pos(vph, 1)

    vph_conf = {'wl_range': [8372.0, 8634.0, 8882.0]}
    vph = create_vph_by_data('VPH863_HR',
                             'HR-I',
                             'v03/VPH863_HR.txt',
                             'v03/VPH863_HR_res.dat',
                             'v03/tvph_0.1aa.dat',
                             conf=vph_conf
                             )
    wheel.put_in_pos(vph, 2)
    return wheel


def create_vph_by_data(name, setup, distortion, resolution, transmission, conf):
    trans = EfficiencyFile(transmission)
    res = EfficiencyFile(resolution)
    vph = MegaraVPH(name=name, setup=setup,
                    vphtable=distortion,
                    resolution=res,
                    transmission=trans,
                    conf=conf)

    return vph


def create_optics():
    _logger.info('create internal optics')
    t = EfficiencyFile('v03/tspect_0.1aa.dat')
    i = InternalOptics(transmission=t)
    return i


def illum1(x, y):
    """Explicit illumination in the focal plane"""
    r = np.hypot(x, y)
    return np.where(r <= 50.0, 1.0, 0.0)


def illum2(x, y):
    """Explicit illumination in the focal plane"""
    r = np.hypot(x, y)
    return 1.0 / (1 + np.exp((r - 130.0) / 10.0))


def create_instrument(insconf):
    # eq = np.random.normal(loc=0.80, scale=0.01, size=(4096,4096))
    # eq = np.clip(eq, 0.0, 1.0)
    # fits.writeto('eq.fits', eq, clobber=True)
    # eq = fits.getdata('eq.fits')

    # Assemble instrument

    detector = create_detector()

    # fibers_mos_base = FiberMOS('FiberMOS')

    # print(fibers_mos_base.config_info())
    cover = MegaraCover()
    focal_plane = FocalPlane(cover)
    fibers_lcb, pseudo_slit_lcb = create_lcb(insconf)
    focal_plane.connect_lcb(fibers_lcb)

    fiber_mos, pseudo_slit_mos = create_mos()
    focal_plane.connect_fibermos(fiber_mos)

    pseudo_slit = PseudoSlitSelector(name='PseudoSlit', capacity=2)
    pseudo_slit.put_in_pos(pseudo_slit_lcb, 0)
    pseudo_slit.put_in_pos(pseudo_slit_mos, 1)

    wheel = create_wheel()
    internal = create_optics()
    _logger.info('create instrument')
    instrument = MegaraInstrument(focal_plane=focal_plane,
                                  wheel=wheel,
                                  pseudo_slit=pseudo_slit,
                                  internal_optics=internal,
                                  detector=detector,
                                  shutter=MegaraShutter()
    )

    cover.set_parent(instrument)
    fiber_mos.set_parent(instrument)
    fibers_lcb.set_parent(instrument)
    return instrument


def create_calibration_unit(illum=None):

    cu = calibrationunit.MegaraCalibrationUnit(capacity=7, name='ICM-MEGARA')

    lamp1 = lamps.BlackBodyLamp('FLAT1', 5400 * u.K, illumination=illum, factor=1e-5)
    lamp2 = lamps.FlatLamp('FLAT2', illumination=illum, factor=5.0)
    lamp3 = lamps.ArcLamp('ARC', illumination=illum)
    lamp4 = lamps.BlackBodyLamp('HALO1', 5400 * u.K, illumination=illum, factor=1e-5)
    lamp5 = lamps.FlatLamp('HALO2', illumination=illum, factor=5.0)
    lamp6 = lamps.ArcLamp('ARC1', illumination=illum)

    cu.put_in_pos('EMPTY', 0)
    cu.put_in_pos(lamp1, 1)
    cu.put_in_pos(lamp2, 2)
    cu.put_in_pos(lamp3, 3)
    cu.put_in_pos(lamp4, 4)
    cu.put_in_pos(lamp5, 5)
    cu.put_in_pos(lamp6, 6)

    return cu


def create_telescope():
    _logger.info('create telescope')
    tel = Telescope(name='GTC', diameter=1040 * u.cm, transmission=EfficiencyFile('v03/ttel_0.1aa.dat'))

    return tel


class ObservingEngine(object):
    def __init__(self, location, zenith_distance):
        self.location = location
        self._zd = zenith_distance

    @property
    def zenith_distance(self):
        return self._zd

    @property
    def airmass(self):
        import math
        za = self._zd
        sec = 1.0 / math.cos(za.to(u.rad).value)
        return sec


def create_observing_engine(location, z):
    _logger.info('create observing engine')
    oe = ObservingEngine(location, z)
    return oe


def restricted_float(x):
    x = float(x)
    if x < 0.0 or x > 36000.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 36000.0]" % (x,))
    return x


if __name__ == '__main__':

    import yaml
    import argparse

    import astropy.coordinates

    from megaradrp.simulation.factory import MegaraImageFactory
    from megaradrp.simulation.atmosphere import AtmosphereModel, ConstSeeing
    from megaradrp.simulation.refraction import DifferentialRefractionModel

    try:
        from numinadb.controldb import ControlSystem1
    except ImportError:
        from megaradrp.simulation.control import ControlSystem

    logging.basicConfig(level=logging.DEBUG)

    parser = argparse.ArgumentParser(prog='megara_sim_b')

    parser.add_argument('-p', '--parameters', metavar="FILE",
                        help="FILE with observing parameters")
    parser.add_argument('-m', '--mos', metavar="FILE",
                        help="FILE with Fiber MOS Configuration")
    parser.add_argument('-t', '--targets', metavar="FILE",
                        help="FILE with target configuration")
    parser.add_argument('-e', '--exposure', type=restricted_float, default=0.0,
                        help="Exposure time per image (in seconds) [0,36000]")
    parser.add_argument('-n', '--nimages', metavar="INT", type=int, default=1,
                        help="Number of images to generate")

    parser.add_argument('omode', choices=megara_sequences().keys(),
                        help="Observing mode of the intrument")

    args = parser.parse_args()

    illum = illum1
    illum = None
    cu = create_calibration_unit(illum=illum)

    import megaradrp.loader as loader
    key = '66f2283e-3049-4d4b-8ef1-14d62fcb611d'
    drp = loader.load_drp()
    insconf = drp.configurations[key]

    instrument = create_instrument(insconf)
    telescope = create_telescope()

    # Observing conditions
    press = 79993.2 * u.Pa
    rel = 0.013333333
    temp = 11.5 * u.deg_C

    # Observing location
    # Madrid coordinates
    location = astropy.coordinates.EarthLocation.from_geodetic(
        lon=-3.7025600,
        lat=40.4165000,
        height=665.0,
        ellipsoid="WGS84"
    )

    refraction_model = DifferentialRefractionModel(temperature=temp, pressure=press, relative_humidity=rel)

    seeing = 0.9
    seeing_model = ConstSeeing(seeing)

    # For twilight spectrum
    factor_tws = 1e-15
    # For night spectrum
    # 4.97490903177e-17 for magnitude 21.9 with  filters/v_johnsonbessel.dat
    factor_ns = 4.97490903177e-17
    atm = AtmosphereModel(twilight=InterpolFile('v03/sky/tw-spec.txt', factor=factor_tws),
                          nightsky=InterpolFile('v03/sky/uves_sky.txt', factor=factor_ns),
                          seeing=seeing_model,
                          extinction=InterpolFile('v03/sky/lapalma_sky_extinction.dat'),
                          refraction=refraction_model
                          )
    telescope.connect(atm)
    factory = MegaraImageFactory()
    observing_engine = create_observing_engine(location, 60 * u.deg)

    control = ControlSystem(factory)
    control.register('MEGARA', instrument)
    control.register('GTC', telescope)
    control.register('ICM-MEGARA', cu)
    control.register('OE', observing_engine)
    control.register('factory', factory)

    # Observation setup
    if args.parameters:
        oparam = yaml.load(open(args.parameters))
        _logger.debug('Configure MEGARA')
        instrument.configure(oparam)
        _logger.debug('Configure ICM-MEGARA')
        cu.configure(oparam)
    if args.mos:
        mosconfig = yaml.load(open(args.mos))
        _logger.debug('Configure Fiber MOS with file %s', args.mos)
        instrument.configure(mosconfig)
    else:
        _logger.debug('Fiber MOS in default positions')

    if args.targets:
        _logger.debug('load targets file %s', args.targets)
        targets = yaml.load(open(args.targets))
    else:
        # empty target list
        targets = dict(central=[0.0, 0.0], targets=dict())

    control.set_targets(targets)

    _logger.info('start simulation')
    control.set_mode(args.omode)

    etime = args.exposure
    repeat = args.nimages

    #yaml.dump(m, open('alldata.yaml', 'rw+'))

    _logger.debug('Exposure time is %f', etime)
    _logger.debug('Number of images is %d', repeat)
    with control:
        control.run(exposure=etime, repeat=repeat)

    import json

    obj = instrument.config_info()

    filename = 'alldata.json'
    with open(filename, 'w') as fd:
        fd.write(json.dumps(obj, sort_keys=True, indent=2,
                            separators=(',', ': ')))
