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
import inspect


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
    def __init__(self, name=None, parent=None):
        self.name = name
        self.parent = None
        self.children = []
        self.set_parent(parent)

    def config_info(self):
        return visit(self)

    def configure(self, meta):
        pass

    def set_parent(self, newparent):
        if self.parent:
            self.parent.children.remove(self)
        self.parent = newparent
        if self.parent:
            self.parent.children.append(self)

    def get_properties(self):
        meta = self.init_config_info()
        for key, prop in inspect.getmembers(self.__class__):
            if isinstance(prop, property):
                try:
                    meta[key] = getattr(self, key).value
                except:
                    meta[key] = getattr(self, key)
        return meta

    def init_config_info(self):
        return dict(name=self.name)

    def end_config_info(self, meta):
        if self.children:
            meta['children'] = [child.name for child in self.children]
        return meta


def visit(node, root='', meta=None):
    sep = '.'
    if meta is None:
        meta = {}

    if node.name is None:
        base = 'unknown'
    else:
        base = node.name

    if root != '':
        node_name = root + sep + base
    else:
        node_name = base

    meta[node_name] = node.get_properties()
    submeta = meta
    for child in node.children:
        visit(child, root=node_name, meta=submeta)
    return meta