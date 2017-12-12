#
# Copyright 2017 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Products of the Megara Pipeline"""

import itertools

import numpy


def to_ds9_reg(obj, ds9reg, rawimage=False, numpix=100, fibid_at=0):
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
    ds9reg.write('#\n# insmode: {0}\n'.format(obj.tags['insmode']))
    ds9reg.write('# vph: {0}\n'.format(obj.tags['vph']))
    ds9reg.write('# uuid: {0}\n'.format(obj.uuid))
    colorbox = ['#ff77ff', '#4444ff']
    itercolors = itertools.cycle(colorbox)
    for fiberdict in obj.contents:
        fibid = fiberdict.fibid
        # boxid = fiberdict.boxid
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
                # Alternate colors
                next_color = next(itercolors)
                ds9reg.write(' # color={0}\n'.format(next_color))
                if fibid_at != 0:
                    if x1 <= fibid_at <= x2:
                        ds9reg.write('text {0} {1} {{{2}}} # color=green '
                                     'font="helvetica 10 bold '
                                     'roman"\n'.format((x1+x2)/2,
                                                       (y1+y2)/2, fibid))
