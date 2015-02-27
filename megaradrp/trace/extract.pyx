
from libc.math cimport floor
import numpy
cimport numpy

cdef int wc_to_pix2(double x):
    return <int>(floor(x + 0.5))

ctypedef fused FType:
    double
    float
    int
    long

ctypedef fused IType:
    int
    long

def extract2(FType[:,:] data, IType[:] xx, double[:] bb1, double[:] bb2, double[:] out):

    cdef size_t size = xx.shape[0]
    cdef size_t i
    cdef int pa, pb
    cdef IType x
    cdef double a,b,w,acc
    
    for i in range(size):
        x = xx[i]
        a = bb1[i]
        b = bb2[i]
        pa, pb = wc_to_pix2(a), wc_to_pix2(b)    
        if pa == pb:
            w = b - a
            out[i] = data[pa, x] * w
        else:
            acc = 0
            w = pa + 0.5 - a
            acc += data[pa, x] * w
            w = b - (pb -0.5)
            acc += data[pb, x] * w
            for c in range(pa + 1, pb):
                acc += data[c,x]
            out[x] = acc
    return out

