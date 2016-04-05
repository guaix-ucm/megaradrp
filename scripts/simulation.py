from __future__ import print_function

import os
import logging
import datetime

import numpy as np
from astropy import units as u

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
from megaradrp.simulation.instrument import PseudoSlit
from megaradrp.simulation.focalplane import FocalPlane
from megaradrp.simulation.detector import ReadParams, MegaraDetectorSat
from megaradrp.simulation.vph import MegaraVPH

from numinadb.dal import Session
from numinadb.model import MyOb, Frame, Base, ObFact

from sqlalchemy import create_engine

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
    wheel = VPHWheel(capacity=3, name='wheel')
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

    cu = calibrationunit.MegaraCalibrationUnit(capacity=7, name='megcalib')

    lamp1 = lamps.BlackBodyLamp('FLAT1', 5400 * u.K, illumination=illum)
    lamp2 = lamps.FlatLamp('FLAT2', photons=7598.34893859, illumination=illum)
    lamp3 = lamps.ArcLamp('ARC', illumination=illum)
    lamp4 = lamps.BlackBodyLamp('HALO1', 5400 * u.K, illumination=illum)
    lamp5 = lamps.FlatLamp('HALO2', photons=7598.34893859, illumination=illum)
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
        self.dbname = 'processing.db'
        uri = 'sqlite:///%s' % self.dbname
        engine = create_engine(uri, echo=False)
        self.conn = None
        self.datadir = 'data'
        Session.configure(bind=engine)

    def register(self, name, element):
        self._elements[name] = element

    def get(self, name):
        return self._elements[name]

    def set_mode(self, mode):
        self.mode = mode

    def run(self, exposure, repeat=1):
        # self.initdb()
        if repeat < 1:
            return

        _logger.info('mode is %s', self.mode)
        try:
            thiss = self.seqs[self.mode]
        except KeyError:
            _logger.error('No sequence for mode %s', self.mode)
            raise

        session = Session()
        now = datetime.datetime.now()
        ob = MyOb(instrument=self.ins, mode=self.mode, start_time=now)
        session.add(ob)

        iterf = thiss.run(self, exposure, repeat)
        count = 1
        for final in iterf:
            _logger.info('image %d of %d', count, repeat)
            name = self.imagecount.runstring()
            fitsfile = factory.create(final, name, self)
            _logger.info('save image %s', name)
            fitsfile.writeto(os.path.join(self.datadir, name), clobber=True)
            # Insert into DB
            newframe = Frame()
            newframe.name = name
            ob.frames.append(newframe)
            count += 1


        # Update completion time of the OB when its finished
        ob.completion_time = datetime.datetime.now()
        # Facts

        from numina.core.pipeline import DrpSystem

        drps = DrpSystem()

        this_drp = drps.query_by_name('MEGARA')

        tagger = None
        for mode in this_drp.modes:
            if mode.key == self.mode:
                tagger = mode.tagger
                break

        if tagger:
            current = os.getcwd()
            os.chdir(self.datadir)
            master_tags = tagger(ob)
            os.chdir(current)
            for k in master_tags:
                fact = ObFact(key=k, value=master_tags[k])
                ob.facts.append(fact)

        session.commit()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.imagecount.__exit__(exc_type, exc_val, exc_tb)

    def initdb(self):

        uri = 'sqlite:///%s' % self.dbname
        engine = create_engine(uri, echo=False)

        Base.metadata.create_all(engine)


class AtmosphereModel(object):

    def __init__(self, twfile):
        self.tw_interp = InterpolFile(twfile)

    def twilight_spectrum(self, wl_in):
        """Twilight spectrum"""
        return 5e4 * self.tw_interp(wl_in)


def restricted_float(x):
    x = float(x)
    if x < 0.0 or x > 36000.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 36000.0]"%(x,))
    return x


if __name__ == '__main__':

    import yaml
    import argparse

    from megaradrp.simulation.factory import MegaraImageFactory
    from megaradrp.simulation.atmosphere import AtmosphereModel

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

    # Observation setup
    # This instruction fixes the number of fibers...

    parser = argparse.ArgumentParser(prog='megara_sim_b')

    parser.add_argument('-p', '--parameters', metavar="FILE",
                        help="FILE with observing parameters")

    parser.add_argument('-e', '--exposure', type=restricted_float, default=0.0,
                        help="Exposure time per image (in seconds) [0,36000]")
    parser.add_argument('-n', '--nimages', metavar="INT", type=int, default=1,
                        help="Number of images to generate")

    parser.add_argument('omode', choices=megara_sequences().keys(),
                        help="Observing mode of the intrument")

    args = parser.parse_args()

    if args.parameters:
        oparam = yaml.load(open(args.parameters))

        _logger.debug('Configure MEGARA with profile %s', oparam['description'])
        instrument.configure(oparam)

        cu.select(oparam['lamp'])

    _logger.info('start simulation')
    control.set_mode(args.omode)

    etime = args.exposure
    repeat = args.nimages

    _logger.debug('Exposure time is %f', etime)
    _logger.debug('Number of images is %d', repeat)
    with control:
        control.run(exposure=etime, repeat=repeat)
