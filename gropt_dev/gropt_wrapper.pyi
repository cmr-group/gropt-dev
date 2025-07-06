from View.MemoryView import __pyx_unpickle_Enum
from __future__ import annotations
import builtins as __builtins__
import numpy as np
import time as time
__all__ = ['GroptParams', 'array_prep', 'get_SAFE', 'np', 'set_verbose', 'time']
class GroptParams:
    @staticmethod
    def __new__(type, *args, **kwargs):
        """
        Create and return a new object.  See help(type) for accurate signature.
        """
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @staticmethod
    def __setstate__(*args, **kwargs):
        ...
    def add_SAFE(self, stim_thresh: float = 1.0, new_first_axis: int = 0, demo_params: bool = True, safe_params: dict = None):
        ...
    def add_bvalue(self, target, tol):
        ...
    def add_gmax(self, gmax):
        ...
    def add_moment(self, order, target):
        ...
    def add_obj_identity(self, weight_mod):
        ...
    def add_smax(self, smax):
        ...
    def diff_init(self, dt: float = 0.0004, TE: float = 0.08, T_90: float = 0.003, T_180: float = 0.005, T_readout: float = 0.016):
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
    def get_out(self):
        ...
    def init(self):
        ...
    def solve(self, min_iter: int = 1, n_iter: int = 2000, gamma_x: float = 1.6, ils_tol: float = 0.001, ils_max_iter: int = 20, ils_min_iter: int = 2, ils_sigma: float = 0.0001):
        ...
    def vec_init_simple(self):
        ...
def __reduce_cython__(self):
    ...
def __setstate_cython__(self, __pyx_state):
    ...
def array_prep(A, dtype, linear = True):
    ...
def get_SAFE(G: np.ndarray, dt: float, true_safe: bool = True, new_first_axis: int = 0, demo_params: bool = True, safe_params: dict = None):
    ...
def set_verbose(level):
    ...
__test__: dict = {}
