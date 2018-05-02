#
# Copyright 2016-2017 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import numpy


class FocalPlane(object):

    def __init__(self, cover):

        self.focalbundle = {}
        self.cover = cover
        self.ddtype = [('fibid', 'i4'),('x', 'f4'), ('y', 'f4'), ('cover', 'f4')]

    def set_cover(self, mode):
        """Cover in the focal plane."""

        self.cover.set_mode(mode)

    def connect_lcb(self, lcb):
        self.focalbundle['LCB'] = lcb

    def connect_fibermos(self, mos):
        self.focalbundle['MOS'] = mos

    def get_visible_fibers(self, name):
        fibid, allpos = self.focalbundle[name].fibers_in_focal_plane()

        p1 = self.cover.visible_fibers(fibid, allpos)

        return numpy.array(p1, dtype=self.ddtype)

    def filter_visible_fibers(self, fibid, allpos):

        p1 = self.cover.visible_fibers(fibid, allpos)

        return numpy.array(p1, dtype=self.ddtype)

    def get_all_fibers(self, name):

        fibid, all_pos = self.focalbundle[name].fibers_in_focal_plane()

        p2 = numpy.empty((len(fibid,)), dtype=self.ddtype)
        p2['fibid'] = fibid
        p2['x'] = all_pos[:,0]
        p2['y'] = all_pos[:,1]
        p2['cover'] = 1.0

        return p2
