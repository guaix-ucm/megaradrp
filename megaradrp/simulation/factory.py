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

import json

import astropy.io.fits as fits


class MegaraImageFactory(object):
    CARDS_P = [
        ('OBSERVAT', 'ORM', 'Name of observatory'),
        ('TELESCOP', 'GTC', 'Telescope id.'),
        ('INSTRUME', 'MEGARA', 'Name of the Instrument'),
        ('ORIGIN', 'SIMULATOR_B', 'FITS file originator'),
        ('OSFILTER', False, 'Sort order filter'),
    ]

    def __init__(self):
        pass

    def bun_fib(self, meta, hdr, id, pos):

        meta_bundle = meta.get('fbundle')

        fibs_id = meta_bundle['fibs_id']
        bunds_id = meta_bundle['bunds_id']
        broken_fibs = meta_bundle['inactive_fibs_id']

        for i in bunds_id:
            key = "BUN%03d_P" % i
            hdr[key] = 0 # priority
            key = "BUN%03d_I" % i
            hdr[key] = "unknown"
            key = "BUN%03d_T" % i
            hdr[key] = "UNASSIGNED"

        # This should be in config_info
        for f, p in zip(id, pos):
            # Coordinates
            key = "FIB%03d_X" % f # X
            hdr[key] = p[0]
            key = "FIB%03d_Y" % f # Y
            hdr[key] = p[1]


        for f, b in zip(fibs_id, bunds_id):
            # Coordinates
            key = "FIB%03d_D" % f # DEC
            hdr[key] = 0.0000
            key = "FIB%03d_R" % f # RA
            hdr[key] = 0.0000
            key = "FIB%03d_O" % f # PA
            hdr[key] = 0.0000

            key = "FIB%03d_A" % f # Active
            if f in broken_fibs:
                hdr[key] = False
            else:
                hdr[key] = True

            key = "FIB%03d_B" % f
            hdr[key] = b
        return hdr

    def create_from_instrument(self, data, name, instrument, mode=''):
        meta = instrument.meta()

        pheader = fits.Header(self.CARDS_P)

        pheader['FILENAME'] = name
        # OBS mode
        pheader['OBSMODE'] = mode

        exptime = meta.get('exposed', 0.0)
        pheader['EXPTIME'] = exptime
        pheader['EXPOSED'] = exptime

        hdu1 = fits.PrimaryHDU(data, header=pheader)
        hdul = fits.HDUList([hdu1])
        return hdul

    def create(self, data, name, control):
        meta = control.get('MEGARA').meta()

        meta_megcalib = control.get('megcalib').config_info()



        mode = control.mode
        pheader = fits.Header(self.CARDS_P)

        pheader['FILENAME'] = name
        # OBS mode
        pheader['OBSMODE'] = mode
        # Detector
        meta_det = meta.get('detector', {})

        exptime = meta_det.get('exposed', 0.0)
        pheader['EXPTIME'] = exptime
        pheader['EXPOSED'] = exptime

        # VPH
        meta_vph = meta.get('vph', {})

        vph_name = meta_vph.get('name', 'unknown')
        pheader['VPH'] = vph_name
        pheader['VPHFWHM1'] = 0.0
        pheader['VPHFWHMC'] = 0.0
        pheader['VPHFWHM2'] = 0.0
        pheader['VPHWL1'] = 0.0
        pheader['VPHWLC'] = 0.0
        pheader['VPHWL2'] = 0.0
        # Focus
        focus = meta.get('focus', 0.0)
        pheader['FOCUS'] = focus

        # Focal plane
        meta_fplane = meta.get('fplane', {})
        cover = meta_fplane.get('cover', 'unknown')
        pheader['COVER'] = cover

        # Instrument mode
        meta_pslit = meta.get('pslit')
        insmode = meta_pslit.get('insmode', 'unknown')
        pheader['INSMODE'] = insmode


        #self.bun_fib(mode, meta, data, pheader)

        if mode in ['arc', 'fiberflat', 'focus']:
            # LAMP keywords

            pheader['LAMP'] = meta_megcalib['label']

        hdu1 = fits.PrimaryHDU(data, header=pheader)

        # Bundles and fibers
        # IN a second extension
        hdu2 = fits.ImageHDU()
        meta_bundle = meta.get('fbundle')
        # FIXME: Ugly hack...
        b = control.get('MEGARA').focal_plane.bundle
        id, pos = b[meta_bundle['name']]

        self.bun_fib(meta, hdu2.header, id, pos)
        hdu2.header['EXTNAME'] = 'FIBERS'

        hdul = fits.HDUList([hdu1, hdu2])
        return hdul


class RunCounter(object):
    """Run number counter"""
    def __init__(self, template, last=1):
        self.template = template
        self.last = last

    def runstring(self):
        """Return the run number and the file name."""
        cfile = self.template % self.last
        self.last += 1
        return cfile


class PersistentRunCounter(RunCounter):
    """Persistent run number counter"""
    def __init__(self, template, last=1, pstore='index.json',):

        last = self.load(pstore, last)

        super(PersistentRunCounter, self).__init__(template, last)

        self.pstore = pstore

    def store(self):
        with open(self.pstore, 'w') as pkl_file:
            json.dump(self.last, pkl_file)

    @staticmethod
    def load(pstore, last):
        file_exists = True

        try:
            with open(pstore, 'rb') as pkl_file:
                last = json.load(pkl_file)
        except IOError:
            file_exists = False

        if not file_exists:
            with open(pstore, 'wb') as pkl_file:
                json.dump(last, pkl_file)

        return last

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.store()


if __name__ == '__main__':

    with PersistentRunCounter('r00%04d') as p:
        for i in range(10):
            print (p.runstring())

    with PersistentRunCounter('r00%04d') as p:
        for i in range(10):
            print (p.runstring())