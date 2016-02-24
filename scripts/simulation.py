from __future__ import print_function

import sys

import logging

import numpy as np
from astropy import units as u
from megaradrp.simulation.instrument import MegaraInstrument
from megaradrp.simulation.factory import MegaraImageFactory
from megaradrp.simulation.actions import simulate_fiber_flat_fits, simulate_focus_fits
from megaradrp.simulation.efficiency import EfficiencyFile
from megaradrp.simulation.instrument import InternalOptics
from megaradrp.simulation import lamps

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
    fibers_lcb = FiberBundle("BUNDLE.LCB", fib_ids, bun_ids, transmission=trans)

    pseudo_slit_lcb = PseudoSlit(name="PSLT.LCB")
    pseudo_slit_lcb.connect_fibers(fib_ids, layouttable[:,2])

    focal_plane.connect_fiber_bundle(fibers_lcb, fib_ids, layouttable[:,0:2])
    return focal_plane, fibers_lcb, pseudo_slit_lcb


def create_mos(focal_plane):
    _logger.info('create mos')
    layouttable = np.loadtxt('v02/MOS_spaxel_centers.dat')
    fib_ids = layouttable[:,4].astype('int').tolist()
    bun_ids = layouttable[:,3].astype('int').tolist()
    trans = EfficiencyFile('v02/tfiber_0.1aa_20m.dat')
    fibers_mos = FiberBundle("BUNDLE.MOS", fib_ids, bun_ids, transmission=trans)

    pseudo_slit_mos = PseudoSlit(name="PSLT.MOS")
    pseudo_slit_mos.connect_fibers(fib_ids, layouttable[:,2])

    focal_plane.connect_fiber_bundle(fibers_mos, fib_ids, layouttable[:,0:2])
    return focal_plane, fibers_mos, pseudo_slit_mos


def create_vph():
    _logger.info('create vphs')
    # vph = create_vph_by_data('VPH405_LR',
    #                          'v02/VPH405_LR2-extra.dat',
    #                          'v02/VPH405_LR_res.dat',
    #                          'v02/tvph_0.1aa.dat'
    #                         )

    vph = create_vph_by_data('VPH926_MR',
                              'v02/VPH926_MR.txt',
                              'v02/VPH926_MR_res.dat',
                              'v02/tvph_0.1aa.dat'
                         )

    # vph = create_vph_by_data('VPH863_HR',
    #                          'v02/VPH863_HR.txt',
    #                          'v02/VPH863_HR_res.dat',
    #                          'v02/tvph_0.1aa.dat'
    #                     )

    return vph


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


if __name__ == '__main__':

    logging.basicConfig(level=logging.DEBUG)


    # eq = np.random.normal(loc=0.80, scale=0.01, size=(4096,4096))
    # eq = np.clip(eq, 0.0, 1.0)
    # fits.writeto('eq.fits', eq, clobber=True)
    # eq = fits.getdata('eq.fits')

    detector = create_detector()

    focal_plane = FocalPlane()
    focal_plane, fibers_lcb, pseudo_slit_lcb = create_lcb(focal_plane)
    focal_plane, fibers_mos, pseudo_slit_mos = create_mos(focal_plane)

    pseudo_slit = dict(lcb=pseudo_slit_lcb, mos=pseudo_slit_mos)
    fibers = dict(lcb=fibers_lcb, mos=fibers_mos)

    vph = create_vph()
    internal = create_optics()
    _logger.info('create instrument')
    instrument = MegaraInstrument(focal_plane=focal_plane,
                                  fibers=fibers,
                                  vph=vph,
                                  pseudo_slit=pseudo_slit,
                                  internal_optics=internal,
                                  detector=detector)

    factory = MegaraImageFactory()

    def illum1(x, y):
        """Explicit illumination in the focal plane"""
        r = np.hypot(x, y)
        return np.where(r <= 50.0, 1.0, 0.5)

    def illum2(x, y):
        """Explicit illumination in the focal plane"""
        r = np.hypot(x, y)
        return 1.0 / (1+np.exp((x-130.0)/ 10.0))

    illum = None

    # Simulated arc spectrum
    wl_in = instrument.vph.wltable_interp()

    # This instruction fixes the number of fibers...
    instrument.set_mode('mos')

    lamp1 = lamps.BlackBodyLamp(5400 * u.K, illumination=illum)
    flat_illum1 = lamp1.illumination_in_focal_plane(instrument, lamp1.flux(wl_in))
    lamp2 = lamps.FlatLamp(photons=7598.34893859, illumination=illum)
    flat_illum2 = lamp2.illumination_in_focal_plane(instrument, lamp2.flux(wl_in))
    lamp3 = lamps.ArcLamp(illumination=illum)
    arc_illum1 = lamp3.illumination_in_focal_plane(instrument, lamp3.flux(wl_in))

    # instrument.set_cover('LEFT')
    #
    # iterf = simulate_fiber_flat_fits(factory, instrument, wltable=wl_in,
    #                                  photons_in=flat_illum1, exposure=20.0, repeat=0)
    # for idx, fitsfile in enumerate(iterf):
    #     fitsfile.writeto('flat-l-%d.fits' % idx, clobber=True)
    #     print('fla2t %d done' % idx)
    #
    # instrument.set_cover('RIGHT')
    # iterf = simulate_fiber_flat_fits(factory, instrument, wltable=wl_in,
    #                                  photons_in=flat_illum1, exposure=20.0, repeat=0)
    # for idx, fitsfile in enumerate(iterf):
    #     fitsfile.writeto('flat-r-%d.fits' % idx, clobber=True)
    #     print('fla2t %d done' % idx)
    #
    # instrument.set_cover('UNSET')
    _logger.info('simulation')

    iterf = simulate_fiber_flat_fits(factory, instrument, wltable=wl_in,
                                     photons_in=flat_illum1, exposure=20.0, repeat=1)
    for idx, fitsfile in enumerate(iterf):
        fitsfile.writeto('flat-u-%d.fits' % idx, clobber=True)
        print('fla2t %d done' % idx)
