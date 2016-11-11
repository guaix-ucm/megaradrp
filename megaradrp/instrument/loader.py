
from __future__ import print_function

import pkgutil
import json
import os

from six import StringIO


class PathLoader(object):
    def __init__(self, inspath, compath):
        self.inspath = inspath
        self.compath = compath

    def build_component_fp(self, key):
        fname = 'component-%s.json' % key
        fcomp = open(os.path.join(self.compath, fname))
        return fcomp

    def build_instrument_fp(self, key):
        fname = 'instrument-%s.json' % key
        fcomp = open(os.path.join(self.inspath, fname))
        return fcomp


class DefaultLoader(object):
    def __init__(self, modpath=None):
        self.modpath = modpath

    def build_component_fp(self, key):
        fname = 'component-%s.json' % key
        return self.build_type_fp(fname)

    def build_instrument_fp(self, key):
        fname = 'instrument-%s.json' % key
        return self.build_type_fp(fname)

    def build_type_fp(self, fname):
        data = pkgutil.get_data(self.modpath, fname)
        fcomp = StringIO(data)
        return fcomp


class Loader(DefaultLoader):
    def __init__(self):
        super(Loader, self).__init__("megaradrp.instrument.configs")


class InstrumentConfiguration(object):

    def __init__(self, instrument):
        self.instrument = instrument
        self.name = 'config'
        self.uuid = '1356a71d-2229-4eac-afce-d0bd11e93cb9'
        self.data_start = 0
        self.data_end = 0
        self.components = {}

    @classmethod
    def from_file(cls, fp, loader=None):

        if loader is None:
            loader = Loader()

        contents = json.load(fp)
        if contents['type'] != 'instrument':
            raise ValueError('type is not instrument')
        instrument = contents['name']
        mm = InstrumentConfiguration(instrument)

        mm.name = contents['description']
        mm.uuid = contents['uuid']
        mm.data_start = 0
        mm.data_end = 0
        for cname, cuuid in contents['components'].items():
            fcomp = loader.build_component_fp(cuuid)
            rr = ComponentConfigurations.from_file(fcomp, loader=loader)
            mm.components[cname] = rr
        return mm

    def get(self, path, **kwds):
        # split key
        vals = path.split('.')
        component = vals[0]
        key = vals[1]
        conf = self.components[component]
        return conf.get(key, **kwds)


class ConfigurationEntry(object):
    def __init__(self, values, depends):
        self.values = values
        self.depends = depends

    def get(self, **kwds):
        result = self.values
        for dep in self.depends:
            key = kwds[dep]
            result = result[key]
        return result

    @classmethod
    def from_file(cls, fp):
        contents = json.load(fp)

        if contents['type'] != 'configuration':
            raise ValueError('type is not configuration')

        key = contents['name']
        confs = contents['configurations']
        val = confs[key]
        mm = ConfigurationEntry(val['values'], val['depends'])
        return mm

class ComponentConfigurations(object):

    def __init__(self):
        self.name = 'config'
        self.uuid = '1356a71d-2229-4eac-afce-d0bd11e93cb9'
        self.data_start = 0
        self.data_end = 0
        self.component = 'component'
        self.configurations = {}

    @classmethod
    def from_file(cls, fp, loader):
        contents = json.load(fp)
        mm = ComponentConfigurations()
        if contents['type'] != 'component':
            raise ValueError('type is not component')
        mm.component = contents['name']
        mm.name = contents['description']
        mm.uuid = contents['uuid']
        mm.data_start = 0
        mm.data_end = 0
        for key, val in contents['configurations'].items():
            if 'uuid' in val:
                # remote component
                fp = loader.build_component_fp(val['uuid'])
                mm.configurations[key] = ConfigurationEntry.from_file(fp)
            else:
                mm.configurations[key] = ConfigurationEntry(val['values'], val['depends'])
        return mm

    def get(self, key, **kwds):
        conf = self.configurations[key]
        return conf.get(**kwds)

def build_instrument_config(uuid, loader=None):

    if loader is None:
        loader = Loader()

    fp = loader.build_instrument_fp(uuid)

    mm = InstrumentConfiguration.from_file(fp, loader=loader)
    return mm



if __name__ == '__main__':

    key = "4fd05b24-2ed9-457b-b563-a3c618bb1d4c"
    mm = build_instrument_config(key, loader=Loader())

    print('detector.scan', mm.get('detector.scan'))
    print('pseudoslit.boxes', mm.get('pseudoslit.boxes', insmode='LCB'))
    print('pseudoslit.boxes_positions', mm.get('pseudoslit.boxes_positions', insmode='LCB', vph='LR-I'))
