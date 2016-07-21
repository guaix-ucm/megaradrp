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

    def bun_fib_lcb(self, meta, hdr):

        sky_bundles_in_lcb = [93, 94, 95, 96, 97, 98, 99, 100]

        extract(hdr, meta, ['MEGARA.LCB', 'nfibers'], 'NFIBERS')
        extract(hdr, meta, ['MEGARA.LCB', 'nbundles'], 'NBUNDLES')

        extract(hdr, meta, ['MEGARA.LCB', 'name'], 'INSMODE')
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
            hdr[key] = fibid
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
        extract(hdr, meta, ['MEGARA.MOS', 'name'], 'INSMODE')


        for i in range(1, nbundles + 1):
            rbpath = 'MEGARA.MOS.RoboticPositioner_%d' % i
            extract(hdr, meta, [rbpath, 'priority'], "BUN%03d_P" % i, default=0)
            extract(hdr, meta, [rbpath, 'something1'], "BUN%03d_I" % i, default="unknown")
            extract(hdr, meta, [rbpath, 'something2'], "BUN%03d_T" % i, default="UNASSIGNED")

            fibs_id = extractm(meta, [rbpath, 'bundle', 'fibs_id'])
            inactive_fibs_id = extractm(meta, [rbpath, 'bundle', 'inactive_fibs_id'])
            pos_fibs = extractm(meta, [rbpath, 'pos'])
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

        pheader = fits.Header(self.CARDS_P)
        pheader['FILENAME'] = name
        pheader['OBSMODE'] = control.mode

        # Seqs
        metacontrol = control.config_info()
        extract(pheader, metacontrol, ['ob_data', 'obsid'], 'OBSID', default=0.0)
        extract(pheader, metacontrol, ['ob_data', 'repeat'], 'NNREP', default=0.0)
        extract(pheader, metacontrol, ['ob_data', 'count'], 'NNSEC', default=0.0)

        instrument = control.get('MEGARA')
        meta = instrument.config_info()
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
        extract(pheader, meta, ['MEGARA.Cover', 'label'], 'cover')
        extract(pheader, meta, ['MEGARA.Cover.Left', 'label'], 'cover1')
        extract(pheader, meta, ['MEGARA.Cover.Right', 'label'], 'cover2')
        extract(pheader, meta, ['MEGARA', 'insmode'], 'insmode', default='unknown')

        #self.bun_fib(mode, meta, data, pheader)

        calibration_unit = control.get('megcalib')
        meta_megcalib = calibration_unit.config_info()
        if control.mode in ['arc', 'fiberflat', 'focus']:
            # LAMP keywords
            extract(pheader, meta_megcalib, ['megcalib', 'label'], 'LAMP')

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


def extract(header, meta, path, key, selector=None, default=None):
    m = meta
    if selector is None:
        selector = lambda x: x
    try:
        for part in path:
            m = m[part]
        header[key] = selector(m)
    except KeyError:
        # Keyword missing
        if default is not None:
            header[key] = default


def extractm(meta, path, selector=None):
    m = meta
    if selector is None:
        selector = lambda x: x
    for part in path:
        m = m[part]
    return selector(m)

if __name__ == '__main__':

    with PersistentRunCounter('r00%04d') as p:
        for i in range(10):
            print (p.runstring())

    with PersistentRunCounter('r00%04d') as p:
        for i in range(10):
            print (p.runstring())
