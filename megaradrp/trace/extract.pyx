
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


cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline int int_min(int a, int b): return a if a <= b else b

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
        pb = int_min(wc_to_pix2(b), data.shape[0])
        pa = int_max(0, wc_to_pix2(a))    
        if pa == pb:
            if pa >= 0 and pa < data.shape[0]:
                w = b - a
                out[x] = data[pa, x] * w
        else:
            acc = 0
            if pa >= 0 and pa < data.shape[0]:
                w = pa + 0.5 - a
                acc += data[pa, x] * w
            
            if pb < data.shape[1] and pb >= 0:
                w = b - (pb -0.5)
                acc += data[pb, x] * w
            for c in range(pa + 1, pb):
                acc += data[c,x]
            out[x] = acc
    return out

