# distutils: language = c++
import time

from libcpp.string cimport string
from libcpp.vector cimport vector

import numpy as np
cimport numpy as cnp
cnp.import_array()

cimport c_gropt

# Prepare numpy array for memory view, with as few copies as possible
def array_prep(A, dtype, linear=True):
    if not A.flags['C_CONTIGUOUS']:
        A = np.ascontiguousarray(A)
    
    A = A.astype(dtype, order='C', copy=False)
    
    if linear:
        A = A.ravel()

    return A 

cdef class GroptParams:
    cdef c_gropt.GroptParams c_gparams

    def __init__(self):
        self.c_gparams = c_gropt.GroptParams() 

    def vec_init_simple(self):
        self.c_gparams.vec_init_simple()

    def diff_init_demo(self):
        self.c_gparams.diff_init_demo()

    def add_gmax(self, gmax):
        self.c_gparams.add_gmax(gmax)

    def add_smax(self, smax):
        self.c_gparams.add_smax(smax)

    def add_moment(self, order, target):
        self.c_gparams.add_moment(order, target)

    def add_bvalue(self, target, tol):
        self.c_gparams.add_bvalue(target, tol)

    def add_obj_identity(self, weight_mod):
        self.c_gparams.add_obj_identity(weight_mod)

    def init(self):
        self.c_gparams.init()

    def solve(self, **kwargs):
        if len(kwargs) == 0:
            # Default solve parameters
            self.c_gparams.solve()
        else:
            # Unpack kwargs and call solve with parameters
            min_iter = kwargs.get('min_iter', 1000)
            n_iter = kwargs.get('n_iter', 2000)
            gamma_x = kwargs.get('gamma_x', 1.6)
            ils_tol = kwargs.get('ils_tol', 1e-3)
            ils_max_iter = kwargs.get('ils_max_iter', 20)
            ils_min_iter = kwargs.get('ils_min_iter', 2)
            ils_sigma = kwargs.get('ils_sigma', 1e-4)

            self.c_gparams.solve(min_iter, n_iter, gamma_x, ils_tol, ils_max_iter, ils_min_iter, ils_sigma)

    def get_out(self):
        cdef double *out
        cdef int out_size
        self.c_gparams.get_output(&out, out_size)
        return np.asarray(<cnp.float64_t[:out_size]> out)

    @property
    def N(self):
        return self.c_gparams.N
    @N.setter
    def N(self, val):
        self.c_gparams.N = val

    @property
    def Naxis(self):
        return self.c_gparams.Naxis
    @Naxis.setter
    def Naxis(self, val):
        self.c_gparams.Naxis = val

    @property
    def dt(self):
        return self.c_gparams.dt
    @dt.setter
    def dt(self, val):
        self.c_gparams.dt = val