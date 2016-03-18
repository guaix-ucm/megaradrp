from __future__ import print_function

import logging

import numpy as np
from astropy import units as u
from megaradrp.simulation.actions import megara_sequences
from megaradrp.simulation.instrument import MegaraInstrument
from megaradrp.simulation.factory import MegaraImageFactory

from megaradrp.simulation.efficiency import EfficiencyFile, InterpolFile
from megaradrp.simulation.instrument import InternalOptics
from megaradrp.simulation.wheel import VPHWheel
from megaradrp.simulation import lamps
from megaradrp.simulation import calibrationunit

from megaradrp.simulation.telescope import Telescope
from megaradrp.simulation.fiberbundle import FiberBundle
from megaradrp.simulation.instrument import PseudoSlit
from megaradrp.simulation.focalplane import FocalPlane
from megaradrp.simulation.detector import ReadParams, MegaraDetectorSat
from megaradrp.simulation.vph import MegaraVPH

_logger = logging.getLogger("simulation")

# create detector from data
def create_detector():
    _logger.info('create detector')
    DSHAPE = (2056 * 2, 2048 * 2)
    PSCAN = 50
    OSCAN = 50
    qe = 1.0 * np.ones(DSHAPE)
    dcurrent = 3.0 / 3600

    readpars1 = ReadParams(gain=1.0, ron=2.0, bias=1000.0)
    readpars2 = ReadParams(gain=1.0, ron=2.0, bias=1005.0)

    qe_wl = EfficiencyFile('v02/tccdbroad_0.1aa.dat')

    detector = MegaraDetectorSat(DSHAPE, OSCAN, PSCAN, qe=qe, qe_wl=qe_wl, dark=dcurrent,
                                   readpars1=readpars1, readpars2=readpars2, bins='11')
    return detector


def create_lcb(focal_plane):
    _logger.info('create lcb')
    layouttable = np.loadtxt('v02/LCB_spaxel_centers.dat')
    fib_ids = layouttable[:,4].astype('int').tolist()
    bun_ids = layouttable[:,3].astype('int').tolist()

    trans = EfficiencyFile('v02/tfiber_0.1aa_20m.dat')
    fibers_lcb = FiberBundle("BUNDLE.LCB", fib_ids, bun_ids, static=True, transmission=trans, inactive=[1, 3])

    pseudo_slit_lcb = PseudoSlit(name="PSLT.LCB", insmode='LCB')
    pseudo_slit_lcb.connect_fibers(fib_ids, layouttable[:,2])

    focal_plane.connect_fiber_bundle(fibers_lcb, fib_ids, layouttable[:,0:2])
    return focal_plane, fibers_lcb, pseudo_slit_lcb


def create_mos(focal_plane):
    _logger.info('create mos')
    layouttable = np.loadtxt('v02/MOS_spaxel_centers.dat')
    fib_ids = layouttable[:,4].astype('int').tolist()
    bun_ids = layouttable[:,3].astype('int').tolist()
    trans = EfficiencyFile('v02/tfiber_0.1aa_20m.dat')
    fibers_mos = FiberBundle("BUNDLE.MOS", fib_ids, bun_ids, static=False, transmission=trans, inactive=[1, 3])

    pseudo_slit_mos = PseudoSlit(name="PSLT.MOS", insmode='MOS')
    pseudo_slit_mos.connect_fibers(fib_ids, layouttable[:,2])

    focal_plane.connect_fiber_bundle(fibers_mos, fib_ids, layouttable[:,0:2])
    return focal_plane, fibers_mos, pseudo_slit_mos


def create_wheel():
    _logger.info('create wheel')
    wheel = VPHWheel(capacity=3)
    _logger.info('create vphs')
    vph = create_vph_by_data('VPH405_LR',
                              'v02/VPH405_LR2-extra.dat',
                              'v02/VPH405_LR_res.dat',
                              'v02/tvph_0.1aa.dat'
                             )
    wheel.put_in_pos(vph, 0)
    vph = create_vph_by_data('VPH926_MR',
                              'v02/VPH926_MR.txt',
                              'v02/VPH926_MR_res.dat',
                              'v02/tvph_0.1aa.dat'
                         )
    wheel.put_in_pos(vph, 1)
    vph = create_vph_by_data('VPH863_HR',
                             'v02/VPH863_HR.txt',
                             'v02/VPH863_HR_res.dat',
                              'v02/tvph_0.1aa.dat'
                         )
    wheel.put_in_pos(vph, 2)
    return wheel


def create_vph_by_data(name, distortion, resolution, transmission):
    trans = EfficiencyFile(transmission)
    res = EfficiencyFile(resolution)
    vph = MegaraVPH(name=name, vphtable=distortion,
                    resolution=res,
                    transmission=trans)

    return vph


def create_optics():
    _logger.info('create internal optics')
    t = EfficiencyFile('v02/tspect_0.1aa.dat')
    i = InternalOptics(transmission=t)
    return i


def illum1(x, y):
    """Explicit illumination in the focal plane"""
    r = np.hypot(x, y)
    return np.where(r <= 50.0, 1.0, 0.5)


def illum2(x, y):
    """Explicit illumination in the focal plane"""
    r = np.hypot(x, y)
    return 1.0 / (1 + np.exp((r - 130.0) / 10.0))


def create_instrument():
 # eq = np.random.normal(loc=0.80, scale=0.01, size=(4096,4096))
    # eq = np.clip(eq, 0.0, 1.0)
    # fits.writeto('eq.fits', eq, clobber=True)
    # eq = fits.getdata('eq.fits')

    # Assemble instrument

    detector = create_detector()

    focal_plane = FocalPlane()
    focal_plane, fibers_lcb, pseudo_slit_lcb = create_lcb(focal_plane)
    focal_plane, fibers_mos, pseudo_slit_mos = create_mos(focal_plane)

    pseudo_slit = dict(lcb=pseudo_slit_lcb, mos=pseudo_slit_mos)
    fibers = dict(lcb=fibers_lcb, mos=fibers_mos)

    wheel = create_wheel()
    internal = create_optics()
    _logger.info('create instrument')
    instrument = MegaraInstrument(focal_plane=focal_plane,
                                  fibers=fibers,
                                  wheel=wheel,
                                  pseudo_slit=pseudo_slit,
                                  internal_optics=internal,
                                  detector=detector)
    return instrument


def create_calibration_unit(illum=None):

    cu = calibrationunit.MegaraCalibrationUnit(capacity=4, name='megcalib')

    lamp1 = lamps.BlackBodyLamp('FLAT1', 5400 * u.K, illumination=illum)
    lamp2 = lamps.FlatLamp('FLAT2', photons=7598.34893859, illumination=illum)
    lamp3 = lamps.ArcLamp('ARC', illumination=illum)

    cu.put_in_pos('EMPTY', 0)
    cu.put_in_pos(lamp1, 1)
    cu.put_in_pos(lamp2, 2)
    cu.put_in_pos(lamp3, 3)

    return cu


def create_telescope():

    tel = Telescope(name='GTC', diameter=100.0, transmission=EfficiencyFile('v02/ttel_0.1aa.dat'))

    return tel


class ControlSystem(object):
    """Top level"""
    def __init__(self):
        self._elements = {}
        from megaradrp.simulation.factory import PersistentRunCounter
        self.imagecount = PersistentRunCounter('r00%04d.fits')
        self.mode = 'null'
        self.ins = 'MEGARA'
        self.seqs = megara_sequences()

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
        count = 1
        for final in iterf:
            _logger.info('image %d of %d', count, repeat)
            name = self.imagecount.runstring()
            fitsfile = factory.create(final, name, self)
            _logger.info('save image %s', name)
            fitsfile.writeto(name, clobber=True)
            count += 1

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.imagecount.__exit__(exc_type, exc_val, exc_tb)


class AtmosphereModel(object):

    def __init__(self, twfile):
        self.tw_interp = InterpolFile(twfile)

    def twilight_spectrum(self, wl_in):
        """Twilight spectrum"""
        return 5e4 * self.tw_interp(wl_in)


if __name__ == '__main__':

    logging.basicConfig(level=logging.DEBUG)

    illum = None
    cu = create_calibration_unit(illum=None)
    instrument = create_instrument()
    telescope = create_telescope()
    atm = AtmosphereModel(twfile='v02/tw-spec.txt')
    telescope.connect(atm)
    factory = MegaraImageFactory()

    control = ControlSystem()
    control.register('MEGARA', instrument)
    control.register('GTC', telescope)
    control.register('megcalib', cu)
    control.register('factory', factory)

    # Obervation setup
    # This instruction fixes the number of fibers...
    instrument.set_mode('LCB')
    # set VPH
    instrument.wheel.select('VPH405_LR')
    #cu.select('FLAT1')

    _logger.info('start simulation')
    control.set_mode('twilightflat')

    with control:
        control.run(exposure=20.0, repeat=1)
