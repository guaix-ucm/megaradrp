#
# Copyright 2016-2017 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#


class FiberBundle(object):
    def __init__(self, name, bid, static=True):
        super(FiberBundle, self).__init__()
        # Geometry of the fibers
        self.name = name
        self.bunds_id = bid
        self.static = static
        self.lf = []

    def add_light_fiber(self, lf):
        self.lf.append(lf)

    @property
    def fibs_id(self):
        return [ch.fibid for ch in self.lf]

    @property
    def inactive_fibs_id(self):
        return [ch.fibid for ch in self.lf if ch.inactive]

    @property
    def active_fibs_id(self):
        return [ch.fibid for ch in self.lf if not ch.inactive]

    def config_info(self):
        return {'name': self.name,
                'nfibers': len(self.lf),
                'fibs_id': self.fibs_id,
                'id': self.bunds_id,
                'static': self.static,
                'inactive_fibs_id': self.inactive_fibs_id
                }