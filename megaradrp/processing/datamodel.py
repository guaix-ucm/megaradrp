#
# Copyright 2016 Universidad Complutense de Madrid
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

"""Data model for MEGARA"""

from __future__ import division

import re
import pkgutil

import astropy.io.fits as fits
from six import StringIO
from numina.flow.datamodel import DataModel


class MegaraDataModel(DataModel):
    """Data model of MEGARA.

    Empty for the moment"""

    def __init__(self):
        # Keys
        self._meta = {
            'texp': ('EXPTIME', None),
            'vph': ('VPH', 'undefined'),
            'obsmode': ('OBSMODE', 'undefined'),
            'tstamp': ('TSTAMP', 'undefined')
        }

    def get_imgid(self, img):
        hdr = self.get_header(img)
        if 'UUID' in hdr:
            return 'uuid:{}'.format(hdr['UUID'])
        elif 'DATE-OBS' in hdr:
            return 'dateobs:{}'.format(hdr['DATE-OBS'])
        else:
            return super(MegaraDataModel, self).get_imgid(img)

    def get_fiberconf(self, img):
        # Obtain FIBER extension
        if 'FIBERS' in img:
            # We have a 'fibers' extension
            # Information os there
            hdr_fiber = img['FIBERS'].header
            return read_fibers_extension(hdr_fiber)
        else:
            # Invalid in general
            insmode = img[0].header.get('INSMODE')
            if insmode == 'LCB':
                data = pkgutil.get_data('megaradrp', 'lcb_default_header.txt')

                default_hdr = StringIO(data)
                hdr_fiber = fits.header.Header.fromfile(default_hdr)
                return read_fibers_extension(hdr_fiber)
            elif insmode == 'MOS':
                raise ValueError('No fibers, MOS')
            else:
                # Read fiber info from headers
                raise ValueError('No fibers')

    def gather_info_dframe(self, img):
        with img.open() as hdulist:
            info = self.gather_info_hdu(hdulist)
        return info

    def gather_info_hdu(self, hdulist):
        values = {}
        values['n_ext'] = len(hdulist)
        extnames = [hdu.header.get('extname', '') for hdu in hdulist[1:]]
        values['name_ext'] = ['PRIMARY'] + extnames
        for key, val in self._meta.items():
            values[key] = hdulist[0].header.get(val[0], val[1])

        return values


class FibersConf(object):

    def __init__(self):
        self.name = ""
        self.conf_id = 1
        self.nbundles = 0
        self.nfibers = 0
        self.bundles = {}
        self.fibers = {}

    def sky_fibers(self, valid_only=False):
        result = []
        for bundle in self.bundles.values():
            if bundle.target_type == 'SKY':
                if valid_only:
                    for fib in bundle.fibers.values():
                        if fib.valid:
                            result.append(fib.fibid)
                else:
                    result.extend(bundle.fibers.keys())
        return result

    def conected_fibers(self, valid_only=False):

        if self.name == 'MOS':
            raise ValueError('not working for MOS')

        result = []
        for bundle in self.bundles.values():
            if bundle.target_type != 'SKY':
                if valid_only:
                    for fib in bundle.fibers.values():
                        if fib.valid:
                            result.append(fib)
                else:
                    result.extend(bundle.fibers.values())
        return result

    def inactive_fibers(self):
        result = []
        for fiber in self.fibers.values():
            if fiber.inactive:
                result.append(fiber.fibid)
        return result

    def active_fibers(self):
        result = []
        for fiber in self.fibers.values():
            if not fiber.inactive:
                result.append(fiber.fibid)
        return result

    def valid_fibers(self):
        result = []
        for fiber in self.fibers.values():
            if fiber.valid:
                result.append(fiber.fibid)
        return result

    def invalid_fibers(self):
        result = []
        for fiber in self.fibers.values():
            if not fiber.valid:
                result.append(fiber.fibid)
        return result


class BundleConf(object):
    def __init__(self):
        self.id = 0
        self.target_type = 'UNASSIGNED'
        self.target_priority = 0
        self.target_name = 'unknown'
        self.x_fix = 0
        self.y_fix = 0
        self.pa_fix = 0
        self.x = 0
        self.y = 0
        self.pa = 0


class FiberConf(object):
    def __init__(self):
        self.fibid = 0
        self.inactive = False
        self.valid = True


def read_fibers_extension(hdr):
    conf = FibersConf()
    conf.name = hdr.get('INSMODE', 'LCB')
    conf.conf_id = hdr.get('CONFID', 1)
    conf.nbundles = hdr.get('NBUNDLES', 89)
    conf.nfibers = hdr.get('NFIBERS', 623)
    # Read bundles

    bun_ids = []
    fib_ids = []
    bundles = conf.bundles
    fibers = conf.fibers

    # loop over everything, count BUN%03d_P and FIB%03d_B
    pattern1 = re.compile(r"BUN(\d+)_P")
    pattern2 = re.compile(r"FIB(\d+)_B")
    for key in hdr:
        bun_re = pattern1.match(key)
        fib_re = pattern2.match(key)
        if bun_re:
            bun_idx = int(bun_re.group(1))
            bun_ids.append(bun_idx)
        elif fib_re:
            fib_idx = int(fib_re.group(1))
            fib_ids.append(fib_idx)

    for i in bun_ids:
        bb = BundleConf()
        bb.id = i
        bb.target_priority = hdr["BUN%03d_P" % i]
        bb.target_name = hdr["BUN%03d_I" % i]
        bb.target_type = hdr["BUN%03d_T" % i]
        bb.fibers = {}
        bundles[i] = bb

    for fibid in fib_ids:
        ff = FiberConf()
        ff.fibid = fibid

        # Coordinates
        ff.d = hdr["FIB%03d_D" % fibid]
        ff.r = hdr["FIB%03d_R" % fibid]
        ff.o = hdr["FIB%03d_O" % fibid]
        # Active
        ff.inactive = not hdr["FIB%03d_A" % fibid]

        # Coordinates XY
        ff.x = hdr["FIB%03d_X" % fibid]
        ff.y = hdr["FIB%03d_Y" % fibid]

        ff.b = hdr["FIB%03d_B" % fibid]

        # Validity
        if ff.inactive:
            ff.valid = False
        else:
            ff.valid = hdr.get("FIB%03d_V" % fibid, True)

        bundles[ff.b].fibers[ff.fibid] = ff
        fibers[ff.fibid] = ff

    return conf