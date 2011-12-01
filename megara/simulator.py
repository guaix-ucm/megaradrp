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
class Wheel(object):
    def __init__(self, size):
        self.pos = 0
        self.size = size

        self.elements = [Holder()] * size
        self.meta = TreeDict()

        for pos in range(self.size):
            self.load(pos, Holder())

        self.position(0)

    def load(self, pos, element):
        if not isinstance(element, Holder):
            h = Holder(element)
        else:
            h = element
        self.elements[pos % self.size] = h

    def position(self, pos):
        self.pos = pos % self.size
        self.meta['active.pos'] = self.pos

        this_elem = self.elements[self.pos]

        self.meta['active.name'] = this_elem.name
        self.meta['active.uid'] = this_elem.uid
        
    def current(self):
        return self.elements[self.pos]

class Holder(object):
    def __init__(self, element=None):
        self.name = 'None' if element is None else element.name
        self.uid = 0 if element is None else element.uid
        self.element = element

class Grism(object):
    def __init__(self, name, uid):
        self.name = name
        self.uid = uid
    
class MegaraWheel(Wheel):

    def __init__(self):
        Wheel.__init__(self, 13)

class MegaraSpectrograph(object):
    def __init__(self, shutter, detector):
        self.shutter = shutter
        self.wheel = MegaraWheel()

        for a in range(13):
            g = Grism('a%d' % (a + 1), a)
            self.wheel.load(a, g)
        
        self.wheel.position(pos=0)


        self.detector = detector
        self.parent = None

        self.path = [self.shutter, self.wheel, self.detector]

        self.meta = TreeDict()
        self.meta['name'] = 'A'
        self.meta['source'] = 'MOS' # One of LFB, SFB, MOS
        self.meta['grism'] = ''
        self.meta['imagetype'] = ''
        self.meta['detector'] = self.detector.meta
        self.meta['wheel'] = self.wheel.meta
 
    def light_path(self, ls):

        for oe in self.path:
            ls = oe.light_path(ls)
        return ls

    def grism(self, pos):
        self.wheel.position(pos)

    def expose(self, time):
        self.detector.expose(time)

    def readout(self):
        data = self.detector.readout() 
        if self.parent is None:
            return self.meta, data
        else:
            meta = self.parent.meta
            meta['spec'] = self.meta          
            return meta, data


    def imagetype(self, name):
        self.meta['imagetype'] = name

class Instrument(object):
    def __init__(self, spectrograph):
        self.spectrograph = spectrograph
        self.spectrograph1 = spectrograph
        self.spec1 = spectrograph
        self.spectrograph.parent = self

        self.meta = TreeDict()
        self.meta['name'] = 'MEGARA'
        self.meta['focalstation'] = 'FCASS'
        self.meta['spec1'] = self.spectrograph.meta
 
    def light_path(self, ls):

        for oe in self.path:
            ls = oe.light_path(ls)
        return ls

class Megara(Instrument):
    def __init__(self):
        shutter = Shutter()
        amplifiers=[Amplifier(shape=(slice(0, 4096), slice(0,4096)), gain=1, ron=1, wdepth=66000)]
        detector = CCDDetector(shape=(4096, 4096), amplifiers=amplifiers, bias=100, dark=0.3)
        ms1 = MegaraSpectrograph(shutter, detector)
        Instrument.__init__(self, ms1)

class MegaraImageFactory(object):
    def __init__(self):
        sfile = StringIO(get_data('megara', 'primary.txt'))
        self.p_templ = pyfits.Header(txtfile=sfile)
        del sfile

    def create(self, metadata, data):

        hh = self.p_templ.copy()
        for rr in hh.ascardlist():
            rr.value = interpolate(metadata, rr.value)

	prim = pyfits.PrimaryHDU(data=data, header=hh)
	hl = [prim]

	hdulist = pyfits.HDUList(hl)
        return hdulist

