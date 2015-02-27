#
import numpy

from megaradrp.trace._extract import extract2

def superex(data, borders, out=None):
    
    if data.dtype.byteorder != '=':
        data2 = data.byteswap().newbyteorder()
    else:
        data2 = data
    
    if out is None:
        out = numpy.zeros((len(borders), data.shape[1]), dtype='float')

    xx = numpy.arange(data2.shape[1])

    for idx, (b1, b2) in enumerate(borders):
        bb1 = b1(xx)    
        bb2 = b2(xx)
        extract2(data2, xx, bb1, bb2, out[idx])
    return out
