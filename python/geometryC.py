# -*- coding: utf-8 -*-
"""
Module containing some geometry functions implemented in C for speed.

To make the bSpline3 dynamically-linked library, on Linux or Mac run:
g++ -c -Wall -Werror -fpic bSpline3.cpp
g++ -shared -o bSpline3.so bSpline3.o

Rick Fenrich 7/5/16
"""

import ctypes
import numpy as np

_mod = ctypes.cdll.LoadLibrary("bSpline3.so")

#void bSplineGeo3(double *knots, double *coefs, double *x, double *y,
#		 double *dydx, int nx, int k, int c);
# Define a special type for the 'double *' argument
class DoubleArrayType:
    def from_param(self, param):
        typename = type(param).__name__
        if hasattr(self, 'from_' + typename):
            return getattr(self, 'from_' + typename)(param)
        elif isinstance(param, ctypes.Array):
            return param
        else:
            raise TypeError("Can't convert %s" % typename)

    # Cast from array.array objects
    def from_array(self, param):
        if param.typecode != 'd':
            raise TypeError('must be an array of doubles')
        ptr, _ = param.buffer_info()
        return ctypes.cast(ptr, ctypes.POINTER(ctypes.c_double))

    # Cast from lists/tuples
    def from_list(self, param):
        val = ((ctypes.c_double)*len(param))(*param)
        return val

    from_tuple = from_list

    # Cast from a numpy array
    def from_ndarray(self, param):
        return param.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

#void bSplineGeo3(double *knots, double *coefs, double *x, double *y,
#		 double *dydx, int nx, int k, int c);
DoubleArray1 = DoubleArrayType()
DoubleArray2 = DoubleArrayType()
DoubleArray3 = DoubleArrayType()
DoubleArray4 = DoubleArrayType()
DoubleArray5 = DoubleArrayType()
_bSplineGeometry = _mod.bSplineGeo3
_bSplineGeometry.argtypes = (DoubleArray1, DoubleArray2, DoubleArray3, \
                            DoubleArray4, DoubleArray5, ctypes.c_int,    \
                            ctypes.c_int, ctypes.c_int)
_bSplineGeometry.restype = None

def bSplineGeometry(knots,coefs,x):
  coefsVec = np.hstack((coefs[0,:],coefs[1,:]))
  y = np.empty(x.size)
  dydx = np.empty(x.size)
  _bSplineGeometry(knots,coefsVec,x,y,dydx,x.size,knots.size,coefsVec.size/2)
  return (y, dydx)

#coefs = np.array(([0.0000, 0.0000, 0.1500, 0.1700, 
#        0.1900, 0.2124, 0.2269, 0.2734, 0.3218, 0.3218, 0.3230, 0.3343, 0.3474, 
#        0.4392, 0.4828, 0.5673, 0.6700, 0.6700],[0.3255, 0.3255, 0.3255, 
#        0.3255, 0.3255, 0.3238, 0.2981, 0.2817, 0.2787, 0.2787, 0.2787, 0.2797, 
#        0.2807, 0.2936, 0.2978, 0.3049, 0.3048, 0.3048]))
#knots = np.hstack(([np.zeros(4), np.arange(1.,coefs.size/2-3),  \
#          np.ones(4)*(coefs.size/2-3)]))
#x = np.linspace(0.,0.68,100)
#(y, dydx) = bSplineGeometry(knots,coefs,x)
#
#print y
#print dydx
#
#plt.plot(x,y)



