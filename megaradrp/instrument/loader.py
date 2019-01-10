
from __future__ import print_function

from numina.core.pipelineload import DefaultLoader


class Loader(DefaultLoader):
    """Instrument configuration loader for MEGARA"""
    def __init__(self):
        super(Loader, self).__init__("megaradrp.instrument.configs")
