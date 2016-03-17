#
# Copyright 2016 Universidad Complutense de Madrid
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


import traceback


class Signal(object):
    """Signal used for callbacks."""
    def __init__(self):
        self.callbacks = []

    def connect(self, callback):
        self.callbacks.append(callback)
        return len(self.callbacks) - 1

    def delete(self, idx):
        self.callbacks.pop(idx)

    def emit(self, *args, **kwds):
        for c in self.callbacks:
            try:
                res = c(*args, **kwds)
                # we can use the result value
                # to disable this callback...
                # not yet implemented
            except TypeError:
                traceback.print_exc()


class HWDevice(object):
    def __init__(self, name=None):
        self.name = name

    def config_info(self):
        return {'name': self.name}

    def configure(self, meta):
        pass

    def set_parent(self, newparent):
        if self.parent:
            self.parent.children.remove(self)
        self.parent = newparent
        if self.parent:
            self.parent.children.append(self)
