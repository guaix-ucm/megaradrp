#
# Copyright 2011-2014 Universidad Complutense de Madrid
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

'''Megara simulator'''

from datetime import datetime
from pkgutil import get_data
from StringIO import StringIO

from astropy.io import fits 
import numpy
from numpy.random import normal, poisson
from numina.treedict import TreeDict
from numina.instrument.detector import DAS, CCDDetector, Amplifier
from numina.instrument import Shutter
from numina.instrument.template import interpolate
from numina.astrotime import datetime_to_mjd

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
        self.detector.connect(self.shutter)
        self.das = DAS(detector)
        self.parent = None

        self.path = [self.shutter, self.wheel, self.detector]

        self.meta = TreeDict()
        self.meta['name'] = 'A'
        self.meta['source'] = 'MOS' # One of LFB, SFB, MOS
        self.meta['grism'] = ''
        self.meta['imagetype'] = ''
        self.meta['detector'] = self.detector.meta
        self.meta['das'] = self.das.meta
        self.meta['wheel'] = self.wheel.meta
 
    def grism(self, pos):
        self.wheel.position(pos)

    def acquire(self, time):
        data = self.das.acquire(time)

        if self.parent is None:
            return self.meta, data
        else:
            meta = self.parent.meta
            meta['spec'] = self.meta
            return meta, data

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

    def connect(self, source):
        self.source = source

class MegaraDetector(CCDDetector):
    def __init__(self):
        amplifiers=[Amplifier(shape=(slice(0, 2048), slice(0,4096)), gain=2.4, ron=6, wdepth=66000),
                    Amplifier(shape=(slice(2048, 4096), slice(0,4096)), gain=2.4, ron=6, wdepth=66000)]
        CCDDetector.__init__(self, shape=(4096, 4096), amplifiers=amplifiers, bias=100, dark=0.3)

    def mode(self, name):
        if name.lower() not in ['slow', 'fast', 'engineering']:
            raise ValueError('%s is not a valid speed mode' % name)
        # FIXME: change here the values of Noise and Gain
        # Slow 40s 3e- 1.2 e-/ADU
        # Fast 23s 6e- 2.4 e-/ADU
        # Engineering 12s - 5e-/ADU
        CCDDetector.mode(self, name)

 
class Megara(Instrument):
    def __init__(self):
        shutter = Shutter()
        detector = MegaraDetector()
        ms1 = MegaraSpectrograph(shutter, detector)
        Instrument.__init__(self, ms1)

class MegaraImageFactory(object):
    def __init__(self):
        #sfile = StringIO(get_data('megara', 'primary.txt'))
        with open('/home/spr/devel/megara/megara/primary.txt') as sfile:
            self.p_templ = fits.Header(txtfile=sfile)
        del sfile

    def create(self, metadata, data):

        print 'create factory'
        hh = self.p_templ.copy()
        for rr in hh.ascardlist():
            rr.value = interpolate(metadata, rr.value)

	prim = fits.PrimaryHDU(data=data, header=hh)
	hl = [prim]

	hdulist = fits.HDUList(hl)
        return hdulist

