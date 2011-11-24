#
# Copyright 2011 Sergio Pascual
# 
# This file is part of Megara DRP
# 
# Pontifex is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Pontifex is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Pontifex.  If not, see <http://www.gnu.org/licenses/>.
#

'''Megara simulator'''

from datetime import datetime
from pkgutil import get_data
from StringIO import StringIO

import pyfits
import numpy
from numpy.random import normal, poisson
from numina.treedict import TreeDict
from numina.instrument import CCDDetector, Amplifier
from numina.instrument.template import interpolate

from pontifex.astrotime import datetime_to_mjd

class OpticalElement(object):
    def __init__(self):
        object.__init__(self)

    def light_path(self, ls):
        return ls

class Shutter(OpticalElement):
    def __init__(self):
        OpticalElement.__init__(self)
        self.opened = True

    def open(self):
        self.opened = True

    def close(self):
        self.opened = False

    def light_path(self, ls):
        if self.opened:
            return ls
        else:
            return None

class Instrument(object):
    def __init__(self, shutter, detector):
        self.shutter = shutter
        self.detector = detector

        self.path = [self.shutter, self.detector]

        self.meta = TreeDict()
        self.meta['name'] = 'MEGARA'
        self.meta['focalstation'] = 'FCASS'
        self.meta['grism0'] = 'A'
        self.meta['imagetype'] = ''
        self.meta['detector'] = self.detector.meta
 
    def light_path(self, ls):

        for oe in self.path:
            ls = oe.light_path(ls)
        return ls

    def grism(self, name):
        self.meta['grism0'] = name

    def expose(self, time):
        self.detector.expose(time)

    def readout(self):
        data = self.detector.readout() 
        return self.meta, data

    def imagetype(self, name):
        self.meta['imagetype'] = name

class Megara(Instrument):
    def __init__(self):
        shutter = Shutter()
        amplifiers=[Amplifier(shape=(slice(0, 4096), slice(0,4096)), gain=1, ron=1, wdepth=66000)]
        detector = CCDDetector(shape=(4096, 4096), amplifiers=amplifiers, bias=100, dark=0.3)
        Instrument.__init__(self, shutter, detector)

class MegaraImageFactory(object):
    def __init__(self):
        sfile = StringIO(get_data('megara', 'primary.txt'))
        self.p_templ = pyfits.Header(txtfile=sfile)
        sfile = StringIO(get_data('megara', 'spec1.txt'))
        self.s1_templ = pyfits.Header(txtfile=sfile)
        sfile = StringIO(get_data('megara', 'spec2.txt'))
        self.s2_templ = pyfits.Header(txtfile=sfile)
        del sfile

    def create(self, metadata, data):

        hh = self.p_templ.copy()
        for rr in hh.ascardlist():
            rr.value = interpolate(metadata, rr.value)

	prim = pyfits.PrimaryHDU(header=hh)
	hl = [prim]

        hh = self.s1_templ.copy()
        for rr in hh.ascardlist():
            rr.value = interpolate(metadata, rr.value)
	spec1 = pyfits.ImageHDU(data=data, header=hh, name='Spec1')
        hl.append(spec1)

	hdulist = pyfits.HDUList(hl)
        return hdulist

