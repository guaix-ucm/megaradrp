#
# Copyright 2011-2017 Universidad Complutense de Madrid
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

"""Products of the Megara Pipeline"""

import numpy
import numpy.polynomial.polynomial as nppol

from .structured import BaseStructuredCalibration

class GeometricTrace(object):
    def __init__(self, fibid, boxid, start, stop, fitparms=None):
        self.fibid = fibid
        self.boxid = boxid
        self.start = start
        self.stop = stop
        self.fitparms = fitparms if fitparms is not None else []
        self.polynomial = None
        # Update polynomial
        self._set_polynomial(fitparms)

    @property
    def valid(self):
        if self.fitparms:
            return True
        else:
            return False

    def __getstate__(self):
        state = self.__dict__.copy()
        del state['polynomial']
        return state

    def __setstate__(self, state):
        self.__dict__ = state
        self._set_polynomial(state['fitparms'])

    def _set_polynomial(self, fitparms):
        if fitparms:
            self.polynomial = nppol.Polynomial(self.fitparms)
        else:
            self.polynomial = nppol.Polynomial([0.0])


class TraceMap(BaseStructuredCalibration):
    """Trace map calibration product"""
    def __init__(self, instrument='unknown'):
        super(TraceMap, self).__init__(instrument)
        self.contents = []
        self.boxes_positions = []
        self.global_offset = nppol.Polynomial([0.0])
        self.ref_column = 2000

    def __getstate__(self):
        st = super(TraceMap, self).__getstate__()
        st['contents'] = [t.__getstate__() for t in self.contents]
        st['boxes_positions'] = self.boxes_positions
        st['global_offset'] = self.global_offset.coef
        st['ref_column'] = self.ref_column
        return st

    def __setstate__(self, state):
        super(TraceMap, self).__setstate__(state)
        self.contents = [GeometricTrace(**trace) for trace in state['contents']]
        self.boxes_positions = state.get('boxes_positions', [])
        self.global_offset = nppol.Polynomial(state.get('global_offset', [0.0]))
        self.ref_column = state.get('ref_column', 2000)
        return self

    def to_ds9_reg(self, ds9reg, rawimage=False, numpix=100, fibid_at=0):
        """Transform fiber traces to ds9-region format.

        Parameters
        ----------
        ds9reg : BinaryIO
            Handle to output file name in ds9-region format.
        rawimage : bool
            If True the traces must be generated to be overplotted on
            raw FITS images.
        numpix : int
            Number of abscissae per fiber trace.
        fibid_at : int
            Abscissae where the fibid is shown (default=0 -> not shown).

        """

        # offset between polynomial and image abscissae
        if rawimage:
            ix_offset = 51
        else:
            ix_offset = 1

        # open output file and insert header

        ds9reg.write('# Region file format: DS9 version 4.1\n')
        ds9reg.write(
            'global color=green dashlist=2 4 width=2 font="helvetica 10 '
            'normal roman" select=1 highlite=1 dash=1 fixed=0 edit=1 '
            'move=1 delete=1 include=1 source=1\n')
        ds9reg.write('physical\n')
        ds9reg.write('#\n# insmode: {0}\n'.format(self.tags['insmode']))
        ds9reg.write('# vph: {0}\n'.format(self.tags['vph']))
        ds9reg.write('# uuid: {0}\n'.format(self.uuid))
        colorbox = ['#ff77ff', '#4444ff']
        for fiberdict in self.contents:
            fibid = fiberdict.fibid
            boxid = fiberdict.boxid
            xmin = fiberdict.start
            xmax = fiberdict.stop
            ds9reg.write('#\n# fibid: {0}\n'.format(fibid))
            # skip fibers without trace
            if fiberdict.valid:
                xp = numpy.linspace(start=xmin, stop=xmax, num=numpix)
                yp = fiberdict.polynomial(xp)
                if rawimage:
                    lcut = (yp > 2056.5)
                    yp[lcut] += 100
                for i in range(len(xp) - 1):
                    x1 = xp[i] + ix_offset
                    y1 = yp[i] + 1
                    x2 = xp[i + 1] + ix_offset
                    y2 = yp[i + 1] + 1
                    ds9reg.write('line {0} {1} {2} {3}'.format(x1, y1, x2, y2))
                    ds9reg.write(' # color={0}\n'.format(colorbox[boxid % 2]))
                    if fibid_at != 0:
                        if x1 <= fibid_at <= x2:
                            ds9reg.write('text {0} {1} {{{2}}} # color=green '
                                         'font="helvetica 10 bold '
                                         'roman"\n'.format((x1+x2)/2,
                                                           (y1+y2)/2, fibid))
