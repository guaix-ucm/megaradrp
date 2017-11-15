#
# Copyright 2016-2017 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
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
    def __init__(self, name, parent=None):
        self.name = name
        self.parent = None
        self.children = {}
        self.set_parent(parent)

    def config_info(self):
        return visit(self)

    def set_parent(self, newparent):
        if self.parent:
            del self.parent.children[self.name]
        if newparent:
            if self.name in newparent.children:
                raise TypeError('%s already registered with newparent')
            else:
                newparent.children[self.name] = self
                self.parent = newparent

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
            meta['children'] = self.children.keys()
        return meta

    def get_device(self, path):
        if path == '':
            return None
        keys = path.split('.')
        root = keys[0]
        trail = keys[1:]
        return self.get_device_pair(root, trail)

    def get_device_pair(self, root, trail):
        if root != self.name:
            return None

        if len(trail) == 0:
            return self
        else:
            node = self.children[trail[0]]
            return node.get_device_pair(trail[0], trail[1:])

    def configure_me(self, value):
        for key in value:
            setattr(self, key, value[key])

    def configure(self, info):
        for key, value in info.items():
            node = self.get_device(key)
            if node:
                node.configure_me(value)


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
    for child in node.children.values():
        visit(child, root=node_name, meta=submeta)
    return meta