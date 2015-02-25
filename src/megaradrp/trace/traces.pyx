import cython
cimport cython

import numpy
cimport numpy

from libc.math cimport floor, ceil
from libc.math cimport fabs
#from libc.stdio cimport printf
from libcpp.vector cimport vector

ctypedef fused FType:
    double
    float
    int
    long

cdef extern from "Trace.h" namespace "Numina":
    cdef cppclass Trace:
        Trace() except +
        void push_back(double x, double y) nogil
        double predict(double x) nogil
        vector[double] xtrace
        vector[double] ytrace
        void reverse() nogil

cdef vector[int] local_max(double* mm, size_t n, double background) nogil:
    
    cdef vector[int] result
    cdef size_t i
    
    if mm[0] >= background:
        if mm[0] > mm[1]:
            result.push_back(0)

    for i in range(1, n-1):
        if mm[i] >= background:
            if (mm[i] > mm[i+1]) and (mm[i] > mm[i-1]):
                result.push_back(i)
                
    if mm[n-1] >= background:
        if mm[n-1] > mm[n-2]:
            result.push_back(n-1)
            
    return result

@cython.boundscheck(False)
cdef vector[double] fit_para_equal_spaced(FType[:] dd) nogil:
    
    cdef vector[double] result
    cdef double A, B, C
    C = dd[1]
    B = 0.5 * (dd[2] - dd[0])
    A = 0.5 * (dd[0] + dd[2] - 2 * dd[1])
    
    result.push_back(A)
    result.push_back(B)
    result.push_back(C)
    return result

@cython.boundscheck(False)
cdef vector[double] interp_max_3(FType[:] dd) nogil:
    '''Parabola that passes through 3 points
        
    With X=[-1,0,1]
    '''
    cdef vector[double] result
    cdef vector[double] params
    cdef double A,B,C
    params = fit_para_equal_spaced(dd)
    A = params[0]
    B = params[1]
    C = params[2]
    result.push_back(-B / (2*A))
    result.push_back(C - B * B / (4*A))
    return result

cdef int wc_to_pix(double x) nogil:
    return <int>floor(x + 0.5)

@cython.cdivision(True)
@cython.boundscheck(False)
cdef int colapse_mean(FType[:, :] arr, vector[double]& out) nogil:
    cdef size_t I = arr.shape[0]
    cdef size_t J = arr.shape[1]
    cdef double accum
    cdef size_t i, j
    for i in range(I):
        accum = 0.0
        for j in range(J):
            accum += arr[i, j]
        out[i] = accum / J
    
    return 0

@cython.cdivision(True)
@cython.boundscheck(False)
cdef Trace _internal_tracing(FType[:, :] arr, Trace& trace, double x, double y,
                             size_t step=1, size_t hs=1, 
                             double maxdis=2.0, double background=150.0,
                             int direction=-1) nogil:

    cdef int col = wc_to_pix(x)
    cdef int row = wc_to_pix(y)
    
    cdef size_t pred_pix
    cdef double prediction
    
    cdef int axis = 1
    cdef size_t i

    cdef size_t axis_size = arr.shape[1]
    cdef const char* aa = "aaaaaa"
    
    # Buffer
    cdef size_t regw = 1 + <int>ceil(maxdis)
    cdef size_t buffsize = 2 * regw + 1
    cdef size_t pred_off 
    cdef vector[double] pbuff
    # Peaks
    cdef vector[int] peaks
    cdef double dis, ndis
    cdef size_t ipeak
    cdef size_t nearp    
    cdef vector[double] result
    
    # Init pbuff
    for _ in range(buffsize):
        pbuff.push_back(0.0)
    
    #print 'Buffer size', buffsize
    
    while (col - step > hs) and (col + step + hs < axis_size):
        #print 'we are in column', col
        col += direction * step
#        printf("--------------------col %i\n", col, prediction)
        #print 'we go to column', col
        #print 'we predict the peak will be in coordinate', prediction
        prediction = trace.predict(col)

        pred_pix = wc_to_pix(prediction)
        pred_off = pred_pix-regw 
        #printf("col %i pred %f pix %i\n", col, prediction, pred_pix)
        # extract a region around the expected peak
        # and collapse it
        colapse_mean(arr[pred_pix-regw:pred_pix+regw + 1,col-hs:col+hs+1], pbuff)

        # Find the peaks
        peaks = local_max(&pbuff[0], buffsize, background=background)
        
        # find nearest peak to prediction
        dis = 40000.0 # a large number
        ipeak = -1
        #printf("npeaks %i\n", peaks.size())
        for i in range(peaks.size()):
            ndis = fabs(peaks[i] + pred_off - prediction)
            if ndis < dis:
                dis = ndis
                ipeak = i
        # check the peak is not further than npixels'
        if ipeak < 0 or dis > maxdis:
            # printf("peak is not found %i\n", ipeak)
            # peak is not found
            return trace
        #printf("peak %f\n", peaks[ipeak] + pred_off)
        #printf("col %i pred %f\n", col, prediction)

        nearp = peaks[ipeak] + pred_pix - regw
        

        # fit the peak with three points
        result = interp_max_3(arr[nearp-1:nearp+2, col])
        
        trace.push_back(col, result[0] + nearp)

    return trace

@cython.cdivision(True)
@cython.boundscheck(False)
def tracing(FType[:, :] arr, double x, double y, size_t step=1, 
                     size_t hs=1, double background=150.0,
                     double maxdis=2.0):
    
    cdef Trace trace 
    # Initial values    
    trace.push_back(x, y)

    _internal_tracing(arr, trace, x, y, step=step, hs=hs, 
                      maxdis=maxdis, background=background,
                      direction=-1)
    trace.reverse()
    _internal_tracing(arr, trace, x, y, step=step, hs=hs, 
                     maxdis=maxdis, background=background,
                     direction=+1)

    result = numpy.empty((trace.xtrace.size(), 2), dtype='float')
    
    for i in range(trace.xtrace.size()):
        result[i,0] = trace.xtrace[i]
        result[i,1] = trace.ytrace[i]
    
    return result

