#
# Copyright 2016-2019 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#


from numina.instrument.components.wheel import Carrousel


class PseudoSlit(object):
    def __init__(self, name, insmode):

        super(PseudoSlit, self).__init__()
        self.name = name
        self.insmode = insmode
        self.fiberset = None
        self.positions = []

    def connect_fibers(self, fiberset, positions):
        self.fiberset = fiberset
        self.positions = positions

    def y_pos(self, fibsid):
        result = []
        for fibid in fibsid:
            pos = self.positions[fibid-1]
            result.append(pos)
        return result

    def config_info(self):
        return {'name': self.name,
                'insmode': self.insmode,
                'nfibers': len(self.fiberset.fibers),
                }


class PseudoSlitSelector(Carrousel):
    pass