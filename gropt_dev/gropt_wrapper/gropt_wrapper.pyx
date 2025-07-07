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

    
    def diff_init(self,
                  dt: float = 400e-6,
                  TE: float = 80e-3,
                  T_90: float = 3e-3,
                  T_180: float = 5e-3,
                  T_readout: float = 16e-3):
        """
        Initialize diffusion sequence parameters.

        Parameters
        ----------
        dt : float, optional
            Raster time in seconds.
        TE : float, optional
            Echo time in seconds.
        T_90 : float, optional
            Duration of the excitation RF pulse in seconds.  This should only
            include the time from excitation, so half of the full RF pulse duration.
        T_180 : float, optional
            Duration of the refocusing RF pulse in seconds.
        T_readout : float, optional
            Time to TE of the readout in seconds.
        """
        
        self.c_gparams.diff_init(dt, TE, T_90, T_180, T_readout)

    def add_gmax(self, gmax):
        self.c_gparams.add_gmax(gmax)

    def add_smax(self, smax):
        self.c_gparams.add_smax(smax)

    def add_moment(self, order, target):
        self.c_gparams.add_moment(order, target)

    def add_SAFE(self, 
                 stim_thresh: float = 1.0,
                 new_first_axis: int = 0, 
                 demo_params: bool = True, 
                 safe_params: dict = None
    ):

        if safe_params is None and not demo_params:
            raise ValueError("If safe_params is None, demo_params must be True.")
        
        cdef double[::1] tau1_view
        cdef double[::1] tau2_view
        cdef double[::1] tau3_view
        cdef double[::1] a1_view
        cdef double[::1] a2_view
        cdef double[::1] a3_view
        cdef double[::1] stim_limit_view
        cdef double[::1] g_scale_view

        cdef double *_unused = NULL

        if safe_params is not None:
            tau1_view = array_prep(safe_params['tau1'], np.float64)
            tau2_view = array_prep(safe_params['tau2'], np.float64)
            tau3_view = array_prep(safe_params['tau3'], np.float64)
            a1_view = array_prep(safe_params['a1'], np.float64)
            a2_view = array_prep(safe_params['a2'], np.float64)
            a3_view = array_prep(safe_params['a3'], np.float64)
            stim_limit_view = array_prep(safe_params['stim_limit'], np.float64)
            g_scale_view = array_prep(safe_params['g_scale'], np.float64)

            self.c_gparams.add_SAFE(stim_thresh,
                                    &tau1_view[0], &tau2_view[0], &tau3_view[0],
                                    &a1_view[0], &a2_view[0], &a3_view[0],
                                    &stim_limit_view[0], &g_scale_view[0],
                                    new_first_axis, False)
        else:
            self.c_gparams.add_SAFE(stim_thresh,
                                    _unused, _unused, _unused, _unused, _unused, _unused,
                                    _unused, _unused,
                                     new_first_axis, demo_params)
        
        
        


    def add_bvalue(self, target, tol):
        self.c_gparams.add_bvalue(target, tol)

    def add_TV(self, tv_lam: float = 0.0, weight_in: float = 1.0):
        """
        Add total variation regularization parameters.

        Parameters
        ----------
        tv_lam : float, optional
            Regularization strength for total variation.
        weight_in : float, optional
            Weighting factor for the input data.
        """
        self.c_gparams.add_TV(tv_lam, weight_in)

    def add_obj_identity(self, weight_mod):
        self.c_gparams.add_obj_identity(weight_mod)

    def init(self):
        self.c_gparams.init()

    def solve(self,
              min_iter: int = 1,
              n_iter: int = 2000,
              gamma_x: float = 1.6,
              ils_tol: float = 1e-3,
              ils_max_iter: int = 20,
              ils_min_iter: int = 2,
              ils_sigma: float = 1e-4,
              ils_tik_lam: float = 0.0):

        self.c_gparams.solve(min_iter, n_iter, gamma_x, ils_tol, ils_max_iter, ils_min_iter, ils_sigma, ils_tik_lam)

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

    @property
    def final_good(self):
        return self.c_gparams.final_good
    @final_good.setter
    def final_good(self, val):
        self.c_gparams.final_good = val

def set_verbose(level):
    c_gropt.set_verbose(level)

def get_SAFE(G: np.ndarray, 
             dt: float, 
             true_safe: bool = True, 
             new_first_axis: int = 0, 
             demo_params: bool = True, 
             safe_params: dict = None
):

    if safe_params is None and not demo_params:
        raise ValueError("If safe_params is None, demo_params must be True.")
    
    cdef double[::1] G_view = array_prep(G, np.float64)

    if G.ndim == 1:
        N = G.size
        Naxis = 1
    else:
        N = G.shape[1]
        Naxis = G.shape[0]

    cdef double *out
    cdef int out_size

    cdef double[::1] tau1_view
    cdef double[::1] tau2_view
    cdef double[::1] tau3_view
    cdef double[::1] a1_view
    cdef double[::1] a2_view
    cdef double[::1] a3_view
    cdef double[::1] stim_limit_view
    cdef double[::1] g_scale_view

    cdef double *_unused = NULL

    if safe_params is not None:
        tau1_view = array_prep(safe_params['tau1'], np.float64)
        tau2_view = array_prep(safe_params['tau2'], np.float64)
        tau3_view = array_prep(safe_params['tau3'], np.float64)
        a1_view = array_prep(safe_params['a1'], np.float64)
        a2_view = array_prep(safe_params['a2'], np.float64)
        a3_view = array_prep(safe_params['a3'], np.float64)
        stim_limit_view = array_prep(safe_params['stim_limit'], np.float64)
        g_scale_view = array_prep(safe_params['g_scale'], np.float64)

        c_gropt.get_SAFE(N, Naxis, dt, &G_view[0], true_safe, new_first_axis, False,
                       &tau1_view[0], &tau2_view[0], &tau3_view[0],
                       &a1_view[0], &a2_view[0], &a3_view[0],
                       &stim_limit_view[0], &g_scale_view[0],
                       &out, out_size)
    else:
        c_gropt.get_SAFE(N, Naxis, dt, &G_view[0], true_safe, new_first_axis, demo_params,
                        _unused, _unused, _unused, 
                        _unused, _unused, _unused,
                        _unused, _unused,
                        &out, out_size)

    
    return np.asarray(<cnp.float64_t[:out_size]> out)