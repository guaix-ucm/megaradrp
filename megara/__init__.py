
__version__ = '0.1.0'

from .simulator import MegaraImageFactory as ImageFactory
from .simulator import Megara as Instrument

__all__ = ['Instrument', 'ImageFactory']
