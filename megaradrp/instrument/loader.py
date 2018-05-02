
from __future__ import print_function

from numina.core.pipelineload import DefaultLoader, build_instrument_config


class Loader(DefaultLoader):
    """Instrument configuration loader for MEGARA"""
    def __init__(self):
        super(Loader, self).__init__("megaradrp.instrument.configs")


if __name__ == '__main__':

    key = "4fd05b24-2ed9-457b-b563-a3c618bb1d4c"
    mm = build_instrument_config(key, loader=Loader())

    print('detector.scan', mm.get('detector.scan'))
    print('pseudoslit.boxes', mm.get('pseudoslit.boxes', insmode='LCB'))
    print('pseudoslit.boxes_positions', mm.get('pseudoslit.boxes_positions', insmode='LCB', vph='LR-I'))
