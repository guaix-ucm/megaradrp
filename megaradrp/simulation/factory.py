#
# Copyright 2015-2018 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Simple monocromatic simulation"""

import uuid
import math
from datetime import datetime

import astropy.wcs as wcs
import astropy.io.fits as fits
from numina.instrument.simulation.factory import extract, extractm


class MegaraImageFactory(object):
    CARDS_P = [
        ('OBSERVAT', 'ORM', 'Name of observatory'),
        ('TELESCOP', 'GTC', 'Telescope id.'),
        ('INSTRUME', 'MEGARA', 'Name of the Instrument'),
        ('ORIGIN', 'SIMULATOR_B', 'FITS file originator'),
        ('OSFILTER', False, 'Sort order filter'),
        ('INSCONF', '66f2283e-3049-4d4b-8ef1-14d62fcb611d')
    ]

    def __init__(self):
        pass

    def bun_fib_lcb(self, meta, hdr):

        sky_bundles_in_lcb = [93, 94, 95, 96, 97, 98, 99, 100]

        extract(hdr, meta, ['MEGARA.LCB', 'nfibers'], 'NFIBERS')
        extract(hdr, meta, ['MEGARA.LCB', 'nbundles'], 'NBUNDLES')
        extract(hdr, meta, ['MEGARA.LCB', 'conf_id'], 'CONFID')
        extract(hdr, meta, ['MEGARA.LCB', 'name'], 'INSMODE')

        # insert WCS
        w = wcs.WCS(naxis=2)
        w.wcs.crpix = [0, 0]
        w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
        w.wcs.cdelt = [1.0 / 3600.0, 1.0 / 3600.0]
        w.wcs.crval = [0.0, 70.0]
        # Rotation around (0,0)
        ang = 1.0 / 180 * math.pi
        cs = math.cos(ang)
        ss = math.sin(ang)
        w.wcs.pc = [[cs, -ss], [ss, cs]]
        # Final rotation to a differente center
        # w.wcs.pv = [[], []]
        hdrwcs = w.to_header()
        hdr.extend(hdrwcs.cards)

        fibers_info = extractm(meta, ['MEGARA.LCB', 'fibers'])
        # FIXME: inactive
        inactive_fibs_id = []
        written_bunds = []
        for fiber_info in fibers_info:
            fibid, bunid, x, y, inactive = fiber_info

            if bunid not in written_bunds:
                written_bunds.append(bunid)

                key = "BUN%03d_P" % bunid
                hdr[key] = 0  # priority
                key = "BUN%03d_I" % bunid
                hdr[key] = "unknown"
                key = "BUN%03d_T" % bunid
                if bunid in sky_bundles_in_lcb:
                    hdr[key] = "SKY"
                else:
                    hdr[key] = "UNASSIGNED"

            key = "FIB%03d_B" % fibid
            hdr[key] = bunid
            # Coordinates
            key = "FIB%03d_D" % fibid  # DEC
            hdr[key] = 0.0000
            key = "FIB%03d_R" % fibid  # RA
            hdr[key] = 0.0000
            key = "FIB%03d_O" % fibid  # PA
            hdr[key] = 0.0000

            key = "FIB%03d_A" % fibid  # Active
            if inactive == 1:
                hdr[key] = False
            else:
                hdr[key] = True

            # Coordinates
            key = "FIB%03d_X" % fibid  # X
            hdr[key] = x
            key = "FIB%03d_Y" % fibid  # Y
            hdr[key] = y

        return hdr

    def bun_fib_mos(self, meta, hdr):

        nbundles = extractm(meta, ['MEGARA.MOS', 'nbundles'])

        extract(hdr, meta, ['MEGARA.MOS', 'nfibers'], 'NFIBERS')
        extract(hdr, meta, ['MEGARA.MOS', 'nbundles'], 'NBUNDLES')
        extract(hdr, meta, ['MEGARA.MOS', 'conf_id'], 'CONFID')
        extract(hdr, meta, ['MEGARA.MOS', 'name'], 'INSMODE')

        for i in range(1, nbundles + 1):
            rbpath = 'MEGARA.MOS.RoboticPositioner_%d' % i
            extract(hdr, meta, [rbpath, 'target_priority'], "BUN%03d_P" % i, default=0)
            extract(hdr, meta, [rbpath, 'target_name'], "BUN%03d_I" % i, default="unknown")
            extract(hdr, meta, [rbpath, 'target_type'], "BUN%03d_T" % i, default="UNASSIGNED")

            fibs_id = extractm(meta, [rbpath, 'bundle', 'fibs_id'])
            inactive_fibs_id = extractm(meta, [rbpath, 'bundle', 'inactive_fibs_id'])
            pos_fibs = extractm(meta, [rbpath, 'fibers_pos'])
            for fibid, pos in zip(fibs_id, pos_fibs):
                extract(hdr, meta, [rbpath, 'id'], "FIB%03d_B" % fibid)
                # Coordinates
                key = "FIB%03d_D" % fibid # DEC
                hdr[key] = 0.0000
                key = "FIB%03d_R" % fibid # RA
                hdr[key] = 0.0000
                key = "FIB%03d_O" % fibid # PA
                hdr[key] = 0.0000

                key = "FIB%03d_A" % fibid # Active
                if fibid in inactive_fibs_id:
                    hdr[key] = False
                else:
                    hdr[key] = True

                # Coordinates
                key = "FIB%03d_X" % fibid # X
                hdr[key] = pos[0]
                key = "FIB%03d_Y" % fibid # Y
                hdr[key] = pos[1]

        return hdr

    def create_from_instrument(self, data, name, instrument, mode=''):
        meta = instrument.config_info()
        pheader = fits.Header(self.CARDS_P)

        # pheader['FILENAME'] = name
        # OBS mode
        pheader['OBSMODE'] = mode
        pheader['DATE-OBS'] = datetime.utcnow().isoformat()
        exptime = meta[instrument.name].get('exposed', 1.0)
        pheader['EXPTIME'] = exptime
        pheader['EXPOSED'] = exptime

        hdu1 = fits.PrimaryHDU(data, header=pheader)
        hdul = fits.HDUList([hdu1])
        return hdul

    def create(self, data, name, control):

        pheader = fits.Header(self.CARDS_P)
        # pheader['FILENAME'] = name
        pheader['OBSMODE'] = control.mode
        pheader['UUID'] = str(uuid.uuid4())
        # Date of simulation
        pheader['DATE'] = datetime.utcnow().isoformat()
        # Date of simulated observation, not set yet
        pheader['DATE-OBS'] = datetime.utcnow().isoformat()

        # Seqs
        metacontrol = control.config_info()
        extract(pheader, metacontrol, ['ob_data', 'obsid'], 'OBSID', default=0.0)
        extract(pheader, metacontrol, ['ob_data', 'repeat'], 'NNREP', default=0.0)
        extract(pheader, metacontrol, ['ob_data', 'count'], 'NNSEC', default=0.0)

        instrument = control.get('MEGARA')
        telescope = control.get('GTC')

        meta = instrument.config_info()
        meta.update(telescope.config_info())

        # not yet implemented
        # pheader['AMPLAYOU'] = "NORMAL"
        # pheader['AMPUP'] = "G"
        # pheader['AMPLOW'] = "E"
        # pheader['GAINUP'] = meta_det.get('gainup', 1.0)
        # pheader['GAINLOW'] = meta_det.get('gainlow', 1.0)
        # pheader['RONUP'] = meta_det.get('ronup', 1.0)
        # pheader['RONLOW'] = meta_det.get('ronlow', 1.0)
        extract(pheader, meta, ['MEGARA.Detector', 'exposed'], 'EXPTIME')
        extract(pheader, meta, ['MEGARA.Detector', 'exposed'], 'EXPOSED')
        extract(pheader, meta, ['MEGARA.Detector', 'vbin'], 'VBIN')
        extract(pheader, meta, ['MEGARA.Detector', 'hbin'], 'HBIN')

        extract(pheader, meta, ['MEGARA.Wheel', 'selected', 'setup'], 'VPH', default='unknown')
        pheader['VPHFWHM1'] = 0.0
        pheader['VPHFWHMC'] = 0.0
        pheader['VPHFWHM2'] = 0.0
        extract(pheader, meta, ['MEGARA.Wheel', 'selected', 'wl_range'], 'VPHWL1', selector=lambda x: x[0])
        extract(pheader, meta, ['MEGARA.Wheel', 'selected', 'wl_range'], 'VPHWLC', selector=lambda x: x[1])
        extract(pheader, meta, ['MEGARA.Wheel', 'selected', 'wl_range'], 'VPHWL2', selector=lambda x: x[2])

        extract(pheader, meta, ['MEGARA.Focus', 'focus'], 'FOCUS', default=0)
        extract(pheader, meta, ['GTC.Focus', 'focus'], 'FOCUST', default=3)
        extract(pheader, meta, ['MEGARA.Cover', 'label'], 'cover')
        extract(pheader, meta, ['MEGARA.Cover.Left', 'label'], 'cover1')
        extract(pheader, meta, ['MEGARA.Cover.Right', 'label'], 'cover2')
        extract(pheader, meta, ['MEGARA', 'insmode'], 'insmode', default='unknown')

        #self.bun_fib(mode, meta, data, pheader)

        calibration_unit = control.get('ICM-MEGARA')
        meta_megcalib = calibration_unit.config_info()
        if control.mode in ['arc', 'fiberflat', 'focus']:
            # LAMP keywords
            extract(pheader, meta_megcalib, ['ICM-MEGARA', 'label'], 'LAMP')

        hdu1 = fits.PrimaryHDU(data, header=pheader)

        # Bundles and fibers
        # IN a second extension
        hdu2 = fits.ImageHDU()

        hdu2.header['EXTNAME'] = 'FIBERS'
        if pheader['insmode'] == 'MOS':
            self.bun_fib_mos(meta, hdu2.header)
        else:
            self.bun_fib_lcb(meta, hdu2.header)

        hdul = fits.HDUList([hdu1, hdu2])
        return hdul
