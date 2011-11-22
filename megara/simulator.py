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
from numina.instrument.template import interpolate

from pontifex.astrotime import datetime_to_mjd

class Amplifier(object):
    def __init__(self):
        self.gain = 1.0
        self.ron = 0.3
        self.shape = (slice(0, 100), slice(0, 100))

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

class CCDDetector(OpticalElement):
    def __init__(self):
        OpticalElement.__init__(self)
        self.shape = (100, 100)
        self.amplifiers = [Amplifier()]
        self.bias = 100.0
        self.dark = 0.1
        self.buffer = numpy.zeros(self.shape)
        self.meta = TreeDict()
        self.light_source = None

        self.meta['readmode'] = 'fast'
        self.meta['readscheme'] = 'perline'
        self.meta['exposed'] = 0
        self.meta['gain'] = 3.6
        self.meta['readnoise'] = 2.16

    def reset(self):
        self.buffer.fill(0)

    def light_path(self, ls):
        self.light_source = ls
        return None

    def expose(self, exposure):
        now = datetime.now()
        # Recording time of start of exposure
        self.meta['exposed'] = exposure
        self.meta['dateobs'] = now.isoformat()
        self.meta['mjdobs'] = datetime_to_mjd(now)

        if self.light_source is not None:
            self.buffer += self.light_source * exposure

        self.buffer = poisson(self.buffer)
        self.buffer += self.dark * exposure

    def readout(self):
        data = self.buffer.copy()
        for amp in self.amplifiers:            
            if amp.ron > 0:
                data[amp.shape] = normal(self.buffer[amp.shape], amp.ron)
            data[amp.shape] /= amp.gain
        data += self.bias
        data = data.astype('int32')
        # readout destroys data
        self.buffer.fill(0)
        return data

    def mode(self, name):
        pass

class Instrument(object):
    def __init__(self, shutter, detector):
        self.shutter = shutter
        self.detector = detector

        self.path = [self.shutter, self.shutter, self.detector]

        self.meta = TreeDict()
        self.meta['name'] = 'MEGARA'
        self.meta['focalstation'] = 'FCASS'
        self.meta['filter0'] = 'B'
        self.meta['imagetype'] = ''
        self.meta['detector'] = self.detector.meta
 
    def light_path(self, ls):

        for oe in self.path:
            ls = oe.light_path(ls)
        return ls

    def filter(self, name):
        self.meta['filter0'] = name

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
        detector = CCDDetector()
        Instrument.__init__(self, shutter, detector)

class MegaraImageFactory(object):
    def __init__(self):
        sfile = StringIO(get_data('megara', 'primary.txt'))
        self.htempl = pyfits.Header(txtfile=sfile)

    def create(self, metadata):

        hh = self.htempl.copy()

        for rr in hh.ascardlist():
            rr.value = interpolate(metadata, rr.value)
        return hh

