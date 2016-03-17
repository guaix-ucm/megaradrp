from __future__ import print_function

import sys

import numpy as np
from astropy import units as u
import astropy.io.fits as fits
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


# create detector from data
def create_detector(mu=1,sigma=0.1):

    DSHAPE = (2056 * 2, 2048 * 2)
    PSCAN = 50
    OSCAN = 50
    if sigma:
        qe = np.random.normal(mu, sigma, DSHAPE)
    else:
        qe = 1.0 * np.ones(DSHAPE)

    fits.writeto('flats/qe.fits',qe, clobber=True)

    dcurrent = 3.0 / 3600

    readpars1 = ReadParams(gain=1.0, ron=2.0, bias=1000.0)
    readpars2 = ReadParams(gain=1.0, ron=2.0, bias=1005.0)

    qe_wl = EfficiencyFile('v02/tccdbroad_0.1aa.dat')

    detector = MegaraDetectorSat(DSHAPE, OSCAN, PSCAN, qe=qe, qe_wl=qe_wl, dark=dcurrent,
                                   readpars1=readpars1, readpars2=readpars2, bins='11')
    return detector

def create_lcb(focal_plane):

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
    t = EfficiencyFile('v02/tvph_0.1aa.dat')
    vph = MegaraVPH(name='VPH405_LR', vphtable='v02/VPH405_LR2-extra.dat', transmission=t)
    return vph


def create_optics():

    t = EfficiencyFile('v02/tspect_0.1aa.dat')
    i = InternalOptics(transmission=t)
    return i


if __name__ == '__main__':

    # eq = np.random.normal(loc=0.80, scale=0.01, size=(4096,4096))
    # eq = np.clip(eq, 0.0, 1.0)
    # fits.writeto('eq.fits', eq, clobber=True)
    # eq = fits.getdata('eq.fits')

    detector = create_detector(mu=1,sigma=0.1)

    focal_plane = FocalPlane()
    focal_plane, fibers_lcb, pseudo_slit_lcb = create_lcb(focal_plane)
    focal_plane, fibers_mos, pseudo_slit_mos = create_mos(focal_plane)

    pseudo_slit = dict(lcb=pseudo_slit_lcb, mos=pseudo_slit_mos)
    fibers = dict(lcb=fibers_lcb, mos=fibers_mos)

    vph = create_vph()
    internal = create_optics()

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

    illum = illum1

    # Simulated arc spectrum
    wl_in = instrument.vph.wltable_interp()

    # This instruction fixes the number of fibers...
    instrument.set_mode('lcb') #lcb|mos

    instrument._internal_focus_factor = 6

    lamp1 = lamps.BlackBodyLamp(5400 * u.K, illumination=illum)
    flat_illum1 = lamp1.illumination_in_focal_plane(instrument, lamp1.flux(wl_in))
    lamp2 = lamps.FlatLamp(illumination=illum)
    flat_illum2 = lamp2.illumination_in_focal_plane(instrument, lamp2.flux(wl_in))
    lamp3 = lamps.ArcLamp(illumination=illum)
    arc_illum1 = lamp3.illumination_in_focal_plane(instrument, lamp3.flux(wl_in))

    # instrument.set_cover('LEFT')

    # iterf = simulate_fiber_flat_fits(factory, instrument, wltable=wl_in,
    #                                  photons_in=flat_illum1, exposure=20.0, repeat=0)
    # for idx, fitsfile in enumerate(iterf):
    #     fitsfile.writeto('flat-l-%d.fits' % idx, clobber=True)
    #     print('fla2t %d done' % idx)

    # instrument.set_cover('RIGHT')
    # iterf = simulate_fiber_flat_fits(factory, instrument, wltable=wl_in,
    #                                  photons_in=flat_illum1, exposure=20.0, repeat=0)
    # for idx, fitsfile in enumerate(iterf):
    #     fitsfile.writeto('flat-r-%d.fits' % idx, clobber=True)
    #     print('fla2t %d done' % idx)

    instrument.set_cover('UNSET')

    iterf = simulate_fiber_flat_fits(factory, instrument, wltable=wl_in,
                                     photons_in=flat_illum1, exposure=100.0, repeat=100)
    for idx, fitsfile in enumerate(iterf):
        fitsfile.writeto('flats/flat-l-%d.fits' % idx, clobber=True)
        print('flats/fla2t %d done' % idx)

    # sys.exit()
    #
    # focii = np.arange(119, 128, 0.5)
    #
    # iterf = simulate_focus_fits(factory, instrument, wl_in, arc_illum1, focii, exposure=20, repeat=1)
    #
    # for idx, fitsfile in enumerate(iterf):
    #     fitsfile.writeto('focus-%d.fits' % idx, clobber=True)
    #     print('focus %d done' % idx)
