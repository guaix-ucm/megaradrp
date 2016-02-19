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


import numpy

import astropy.io.fits as fits


class MegaraImageFactory(object):
    CARDS_P = [
        ('OBSERVAT', 'ORM', 'Name of observatory'),
        ('TELESCOP', 'GTC', 'Telescope id.'),
        ('INSTRUME', 'MEGARA', 'Name of the Instrument'),
        ('ORIGIN', 'Simulator', 'FITS file originator'),
    ]

    def create(self, mode, meta, data):
        pheader = fits.Header(self.CARDS_P)

        # Detector
        meta_det = meta.get('detector', {})

        exptime = meta_det.get('exposed', 0.0)
        pheader['EXPTIME'] = exptime
        pheader['EXPOSED'] = exptime

        # VPH
        meta_vph = meta.get('vph', {})

        vph_name = meta_vph.get('name', 'unknown')
        pheader['VPH'] = vph_name

        # Focus
        focus = meta.get('focus', 0.0)
        pheader['FOCUS'] = focus

        # Focal plane
        meta_fplane = meta.get('fplane', {})
        cover = meta_fplane.get('cover', 'unknown')
        pheader['COVER'] = cover

        hdu1 = fits.PrimaryHDU(data, header=pheader)
        hdul = fits.HDUList([hdu1])
        return hdul
