# distutils: language = c++
import time

from libcpp.string cimport string
from libcpp.vector cimport vector

import numpy as np
cimport numpy as cnp
cnp.import_array()

cimport c_gropt

def array_prep(A, dtype, linear=True):
    """
    Prepare a NumPy array for a Cython memory view.

    This function ensures that the input array `A` is C-contiguous and
    has the specified data type, making it suitable for efficient use
    with Cython memory views. It aims to minimize data copies by using
    `np.ascontiguousarray` and `astype(copy=False)` where possible.
    The array can also be flattened.

    Parameters
    ----------
    A : array_like
        The input array to be prepared.
    dtype : numpy.dtype
        The target data type for the array.
    linear : bool, optional
        If True (default), the array is flattened into a 1D array
        using `ravel()`. If False, the array's dimensions are preserved.

    Returns
    -------
    numpy.ndarray
        The prepared NumPy array, which is C-contiguous, has the correct
        `dtype`, and is optionally flattened.

    """
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


    def vec_init_simple(self,
                        N: int = -1,
                        Naxis: int = -1,
                        first_val: float = 0.0,
                        last_val: float = 0.0):
        """
        Initialize the set_vec and inv_vec settings in gparams.

        Parameters
        ----------
        N : int, optional
            Number of points in a *single* axis of the waveform. Negative values
            will use the existing value.
        Naxis : int, optional
            Number of axes in the waveform. Negative values will use the existing value.
        first_val : float, optional
            Fixed value for the first point in the gradient vector. [mT/m]
        last_val : float, optional
            Fixed value for the last point in the gradient vector. [mT/m]
        """
        self.c_gparams.vec_init_simple(N, Naxis, first_val, last_val)


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


    def setvec_X0(self, 
                  X0: np.ndarray,  
                  set_others: bool = True):
        """
        Set the initial waveform guess.

        Parameters
        ----------
        X0 : np.ndarray
            Initial guess for the waveform (warm start). [T/m]
        set_others : bool, optional
            If True, other related vectors (inv_vec, set_vals, fixer) will be updated
            with standard values.
        """
        cdef double[::1] X0_view = array_prep(X0, np.float64)

        if X0.ndim == 1:
            N = X0.size
            Naxis = 1
        else:
            N = X0.shape[1]
            Naxis = X0.shape[0]

        self.c_gparams.setvec_X0(N, Naxis, &X0_view[0], set_others)


    def set_ils_solver(self, ils_method: str = 'CG'):
        """
        Set the indirect solver method.

        Parameters
        ----------
        ils_method : str
            The name of the indirect solver method to use.
            Currently supported methods are 'CG', 'NLCG', and 'BiCGstabl'. (case-sensitive)
        """
        self.c_gparams.set_ils_solver(ils_method.encode('utf-8'))


    def add_gmax(self, 
                 gmax: float = 0.03,
                 rot_variant: bool = True,
                 weight_mod: float = 1.0):
        """
        Adds a constraint for the maximum gradient amplitude.

        Parameters
        ----------
        gmax : float, optional
            The maximum allowed gradient magnitude [T/m].
            Defaults to 0.03.
        rot_variant : bool, optional
            If True, uses the rotationally invariant formulation of the gmax
            constraint. i.e. each gradient axis can hit gmax.  Defaults to True.
        weight_mod : float, optional
            A weighting factor for this specific constraint in the optimization
            problem. Defaults to 1.0.
        """
        self.c_gparams.add_gmax(gmax, rot_variant, weight_mod)


    def add_smax(self, 
                 smax: float = 80.0,
                 rot_variant: bool = True,
                 weight_mod: float = 1.0):
        """
        Adds a constraint for the maximum gradient slew rate.

        Parameters
        ----------
        smax : float, optional
            The maximum allowed gradient slew rate [T/m/s].
            Defaults to 0.03.
        rot_variant : bool, optional
            If True, uses the rotationally invariant formulation of the smax
            constraint. i.e. each gradient axis can hit smax.  Defaults to True.
        weight_mod : float, optional
            A weighting factor for this specific constraint in the optimization
            problem. Defaults to 1.0.
        """
        self.c_gparams.add_smax(smax, rot_variant, weight_mod)

   
    def add_moment(self, 
                   order: int = 0, 
                   target: float = 0.0, 
                   tol: float = 1e-6, 
                   units: str = 'mT*ms/m',
                   axis: int = 0, 
                   start_idx: int = -1, 
                   stop_idx: int = -1, 
                   ref_idx: int = 0,
                   weight_mod: float = 1.0):
        """
        
        Adds a moment constraint to the optimization problem.

        Parameters
        ----------
        order : int, optional
            The order of the moment.
            Defaults to 0.
        target : float, optional
            The target value for the moment constraint.
            Defaults to 0.0.
        tol : float, optional
            The tolerance for satisfying the moment constraint.
            Defaults to 1e-6.
        units : str, optional
            The units of the moment constraint, for `target` and `tolerance`
            Defaults to 'mT*ms/m'.
            Can be 'mT*ms/m', 'T*s/m', 'rad*s/m', or 's/m'
        axis : int, optional
            The axis along which the moment is calculated.
            Defaults to 0.
        start_idx : int, optional
            The starting index (inclusive) for the range over which the
            moment is calculated. A value of -1 indicates the beginning of the
            axis. 
            Defaults to -1.
        stop_idx : int, optional
            The stopping index (exclusive) for the range over which the
            moment is calculated. A value of -1 indicates the end of the
            axis. 
            Defaults to -1.
        ref_idx : int, optional
            The reference index, (i.e. index where t=0 in moment calculations), 
            which often could be the time of excitation.
            Defaults to 0.
        weight_mod : float, optional
            A weighting factor for this specific constraint in the optimization
            problem. Defaults to 1.0.
        """
        self.c_gparams.add_moment(order, target, tol, units.encode('utf-8'), axis, start_idx, stop_idx, ref_idx, weight_mod)

    def add_SAFE(self, 
                 stim_thresh: float = 1.0,
                 new_first_axis: int = 0, 
                 demo_params: bool = True, 
                 safe_params: dict = None,
                 weight_mod: float = 1.0):
        """
        Adds a SAFE constraint to the optimization problem.

        Parameters
        ----------
        stim_thresh : float
            The stimulus threshold for the SAFE constraint.
            Defaults to 1.0.
        new_first_axis : int
            Swaps the first axis of the SAFE parameters used for the constraint. 
            This is useful for a single-axis optimization, where you want to test
            the SAFE parameters for a different axis than the first.
            Defaults to 0.
        demo_params : bool
            Whether to use demo parameters for the SAFE constraint.
            NOTE: If `safe_params` is not None, this parameter is ignored.
            Defaults to True.
        safe_params : dict, optional
            A dictionary of SAFE parameters. See `gropt.readasc`
        weight_mod : float, optional
            A weighting factor for this specific constraint in the optimization
            problem. Defaults to 1.0.
        """
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
                                    new_first_axis, False, weight_mod)
        else:
            self.c_gparams.add_SAFE(stim_thresh,
                                    _unused, _unused, _unused, _unused, _unused, _unused,
                                    _unused, _unused,
                                     new_first_axis, demo_params, weight_mod)

    def add_SAFE2(self, 
                 stim_thresh: float = 1.0,
                 new_first_axis: int = 0, 
                 demo_params: bool = True, 
                 safe_params: dict = None,
                 weight_mod: float = 1.0):
        """
        Adds a SAFE constraint to the optimization problem.

        Parameters
        ----------
        stim_thresh : float
            The stimulus threshold for the SAFE constraint.
            Defaults to 1.0.
        new_first_axis : int
            Swaps the first axis of the SAFE parameters used for the constraint. 
            This is useful for a single-axis optimization, where you want to test
            the SAFE parameters for a different axis than the first.
            Defaults to 0.
        demo_params : bool
            Whether to use demo parameters for the SAFE constraint.
            NOTE: If `safe_params` is not None, this parameter is ignored.
            Defaults to True.
        safe_params : dict, optional
            A dictionary of SAFE parameters. See `gropt.readasc`
        weight_mod : float, optional
            A weighting factor for this specific constraint in the optimization
            problem. Defaults to 1.0.
        """
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

            self.c_gparams.add_SAFE2(stim_thresh,
                                    &tau1_view[0], &tau2_view[0], &tau3_view[0],
                                    &a1_view[0], &a2_view[0], &a3_view[0],
                                    &stim_limit_view[0], &g_scale_view[0],
                                    new_first_axis, False, weight_mod)
        else:
            self.c_gparams.add_SAFE2(stim_thresh,
                                    _unused, _unused, _unused, _unused, _unused, _unused,
                                    _unused, _unused,
                                     new_first_axis, demo_params, weight_mod)

    def add_SAFE_vec(self, 
                     stim_thresh_vec: np.ndarray,
                     new_first_axis: int = 0, 
                     demo_params: bool = True, 
                     safe_params: dict = None,
                     weight_mod: float = 1.0):
        """
        Adds a SAFE constraint to the optimization problem with a vector stimulation limit.


        Parameters
        ----------
        stim_thresh : float
            The stimulus threshold for the SAFE constraint.
            Defaults to 1.0.
        new_first_axis : int
            Swaps the first axis of the SAFE parameters used for the constraint. 
            This is useful for a single-axis optimization, where you want to test
            the SAFE parameters for a different axis than the first.
            Defaults to 0.
        demo_params : bool
            Whether to use demo parameters for the SAFE constraint.
            NOTE: If `safe_params` is not None, this parameter is ignored.
            Defaults to True.
        safe_params : dict, optional
            A dictionary of SAFE parameters. See `gropt.readasc`
        weight_mod : float, optional
            A weighting factor for this specific constraint in the optimization
            problem. Defaults to 1.0.
        """
        if safe_params is None and not demo_params:
            raise ValueError("If safe_params is None, demo_params must be True.")
        
        cdef double[::1] stim_view = array_prep(stim_thresh_vec, np.float64)
        cdef int N_vec = stim_thresh_vec.size

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

            self.c_gparams.add_SAFE_vec(N_vec, &stim_view[0],
                                    &tau1_view[0], &tau2_view[0], &tau3_view[0],
                                    &a1_view[0], &a2_view[0], &a3_view[0],
                                    &stim_limit_view[0], &g_scale_view[0],
                                    new_first_axis, False, weight_mod)
        else:
            self.c_gparams.add_SAFE_vec(N_vec, &stim_view[0],
                                    _unused, _unused, _unused, _unused, _unused, _unused,
                                    _unused, _unused,
                                     new_first_axis, demo_params, weight_mod)
    
    def add_SAFE2_vec(self, 
                     stim_thresh_vec: np.ndarray,
                     new_first_axis: int = 0, 
                     demo_params: bool = True, 
                     safe_params: dict = None,
                     weight_mod: float = 1.0):
        """
        Adds a SAFE constraint to the optimization problem with a vector stimulation limit.


        Parameters
        ----------
        stim_thresh : float
            The stimulus threshold for the SAFE constraint.
            Defaults to 1.0.
        new_first_axis : int
            Swaps the first axis of the SAFE parameters used for the constraint. 
            This is useful for a single-axis optimization, where you want to test
            the SAFE parameters for a different axis than the first.
            Defaults to 0.
        demo_params : bool
            Whether to use demo parameters for the SAFE constraint.
            NOTE: If `safe_params` is not None, this parameter is ignored.
            Defaults to True.
        safe_params : dict, optional
            A dictionary of SAFE parameters. See `gropt.readasc`
        weight_mod : float, optional
            A weighting factor for this specific constraint in the optimization
            problem. Defaults to 1.0.
        """
        if safe_params is None and not demo_params:
            raise ValueError("If safe_params is None, demo_params must be True.")
        
        cdef double[::1] stim_view = array_prep(stim_thresh_vec, np.float64)
        cdef int N_vec = stim_thresh_vec.size

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

            self.c_gparams.add_SAFE2_vec(N_vec, &stim_view[0],
                                    &tau1_view[0], &tau2_view[0], &tau3_view[0],
                                    &a1_view[0], &a2_view[0], &a3_view[0],
                                    &stim_limit_view[0], &g_scale_view[0],
                                    new_first_axis, False, weight_mod)
        else:
            self.c_gparams.add_SAFE2_vec(N_vec, &stim_view[0],
                                    _unused, _unused, _unused, _unused, _unused, _unused,
                                    _unused, _unused,
                                     new_first_axis, demo_params, weight_mod)


    def add_bvalue(self, 
                   target: float = 100.0, 
                   tol: float = 1.0,
                   start_idx0: int = -1, 
                   stop_idx0: int = -1, 
                   weight_mod: float = 1.0,
                   mode: int = 1,
                   max_scale: float = 1.01):
        """
        Adds a b-value constraint to the optimization problem.
        
        Parameters
        ----------
        target : float, optional
            The target b-value for the constraint. 
            Defaults to 100.0.
        tol : float, optional
            The tolerance for the b-value constraint. 
            Keeping this value a little larger (>0.1) might make the optimization faster.
            Defaults to 1.0.
        start_idx0 : int, optional
            The starting index (inclusive) for the range over which the constraint is applied. 
            Defaults to -1 ( = full waveform).
        stop_idx0 : int, optional
            The stopping index (exclusive) for the range over which the constraint is applied. 
            Defaults to -1 ( = full waveform).
        weight_mod : float, optional
            A weighting factor for this specific constraint in the optimization problem. 
            Defaults to 1.0.
        mode : int, optional
            The mode for the b-value constraint. 
            1 = set b-value to target (default)
            2 = set minimum b-value to target
        max_scale : float, optional
            How much to scale the proximal bvalue when mode=3, maximizing "constraint"
            Defaults to 1.01 .
        """
        self.c_gparams.add_bvalue(target, tol, start_idx0, stop_idx0, weight_mod, mode, max_scale)

    def add_TV(self, 
               tv_lam: float = 0.0, 
               weight_mod: float = 1.0):
        """
        Add total variation regularization parameters.

        Unlike other constraints, this is not required to be feasible.
        TODO: Make this property a user choice so it can be required for feasibility.
        Parameters
        ----------
        tv_lam : float, optional
            Regularization strength for total variation.  This wont do
            anything if it is not > 0.0.
        weight_mod : float, optional
            A weighting factor for this specific constraint in the optimization problem. 
            Defaults to 1.0.
        """
        self.c_gparams.add_TV(tv_lam, weight_mod)

    def add_obj_identity(self, 
                         weight_mod: float = 1.0):
        """
        Add an identity objective function to the problem.

        This is typically used for regularization, penalizing the L2 norm
        of the solution vector.

        Parameters
        ----------
        weight_mod : float
            The weighting factor for the identity objective. Larger values
            increase the penalty on the solution's norm.
            Default 1.0
        """
        self.c_gparams.add_obj_identity(weight_mod)


    def prepare(self):
        """
        Force a preparation of the GroptParams object.

        This mostly just allocates vectors in each of the constraints, and 
        sets up the initial optimization variables.  It is automatically done in solve
        if it has not been done yet.
        """
        self.c_gparams.prepare()

    def solve(self,
              min_iter: int = 1,
              max_iter: int = 2000,
              gamma_x: float = 1.6,
              ils_tol: float = 1e-3,
              ils_max_iter: int = 20,
              ils_min_iter: int = 2,
              ils_sigma: float = 1e-4,
              ils_tik_lam: float = 0.0):
        """
        Run the optimization solver.

        Parameters
        ----------
        min_iter : int, optional
            Minimum number of iterations for the main SDMM loop. The solver will
            run for at least this many iterations. Defaults to 1.
        max_iter : int, optional
            Maximum number of iterations for the main SDMM loop. Defaults to 2000.
        gamma_x : float, optional
            Relaxation parameter for the SDMM update step.
            correspond to over-relaxation, which can speed up convergence.
            Defaults to 1.6.
        ils_tol : float, optional
            Relative tolerance for the inner-loop indirect linear solver (e.g., CG).
            The ILS stops when the residual norm is less than `ils_tol` times
            the initial residual norm. Defaults to 1e-3.
        ils_max_iter : int, optional
            Maximum number of iterations for the indirect linear solver per main
            SDMM iteration. Defaults to 20.
        ils_min_iter : int, optional
            Minimum number of iterations for the indirect linear solver.
            Defaults to 2.
        ils_sigma : float, optional
            Regularization parameter for the linear system solved in the inner
            loop, related to the ADMM penalty parameter. Defaults to 1e-4.
        ils_tik_lam : float, optional
            Tikhonov (L2) regularization parameter for the indirect linear
            solver. Adds a penalty on the norm of the solution to improve
            conditioning. Defaults to 0.0.
        """
        self.c_gparams.solve(min_iter, max_iter, gamma_x, ils_tol, ils_max_iter, ils_min_iter, ils_sigma, ils_tik_lam)

    def get_out(self) -> np.ndarray:
        """
        Get the gradient waveform output from the C++ layer.

        Returns
        -------
        np.ndarray
            The output array.
        """
        cdef double *out
        cdef int out_size
        self.c_gparams.get_output(&out, out_size)
        return np.asarray(<cnp.float64_t[:out_size]> out)

    def get_output_bvalue(self) -> float:
        """
        Get the b-value of the current waveform.

        Returns
        -------
        float
            The b-value of the current waveform.
        """
        return self.c_gparams.get_output_bvalue()

    def test_reduce_and_solve(self):
        """
        Test the reduction and solve process.

        This method will call the `test_reduce_and_solve` method in the C++ layer,
        which will reduce the problem and then solve it.
        """
        self.c_gparams.test_reduce_and_solve()

    # TODO: Make these protected and use getter/setters (except final_good)
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

    @property
    def final_n_feval(self):
        return self.c_gparams.final_n_feval
    @final_n_feval.setter
    def final_n_feval(self, val):
        self.c_gparams.final_n_feval = val

cdef class SolverGroptSDMM:
    cdef c_gropt.SolverGroptSDMM c_solver

    def __init__(self):
        self.c_solver = c_gropt.SolverGroptSDMM()

    def set_general_params(self, 
                        min_iter: int = 1, 
                        max_iter: int = 2000, 
                        log_interval: int = 20, 
                        gamma_x: float = 1.6,
                        max_feval: int = 12000):
        """
        Set general parameters for the SDMM solver.

        Parameters
        ----------
        min_iter : int, optional
            Minimum number of iterations for the main SDMM loop. The solver will
            run for at least this many iterations. Defaults to 1.
        max_iter : int, optional
            Maximum number of iterations for the main SDMM loop. Defaults to 2000.
        log_interval : int, optional
            Interval at which to log the progress of the solver. 
            This will only be seen when verbose logging is enabled.
            Defaults to 20.
        gamma_x : float, optional
            Relaxation parameter for the SDMM update step.
            correspond to over-relaxation, which can speed up convergence.
            Defaults to 1.6.
        """
        self.c_solver.set_general_params(min_iter, max_iter, log_interval, gamma_x, max_feval)


    def set_ils_params(self, 
                    ils_tol: float = 1e-3, 
                    ils_max_iter: int = 20, 
                    ils_min_iter: int = 2, 
                    ils_sigma: float = 1e-4, 
                    ils_tik_lam: float = 0.0):
        """ 
        Set parameters for the indirect linear solver used in the SDMM.

        Parameters
        ----------
        ils_tol : float, optional
            Relative tolerance for the inner-loop indirect linear solver (e.g., CG).
            The ILS stops when the residual norm is less than `ils_tol` times
            the initial residual norm. Defaults to 1e-3.
        ils_max_iter : int, optional
            Maximum number of iterations for the indirect linear solver per main
            SDMM iteration. Defaults to 20.
        ils_min_iter : int, optional
            Minimum number of iterations for the indirect linear solver.
            Defaults to 2.
        ils_sigma : float, optional
            Regularization parameter for the linear system solved in the inner
            loop, related to the ADMM penalty parameter. Defaults to 1e-4.
        ils_tik_lam : float, optional
            Tikhonov (L2) regularization parameter for the indirect linear
            solver. Adds a penalty on the norm of the solution to improve
            conditioning. Defaults to 0.0.
        """
        self.c_solver.set_ils_params(ils_tol, ils_max_iter, ils_min_iter, ils_sigma, ils_tik_lam)

    def set_sdmm_params(self, 
                        rw_interval: int = 8, 
                        rw_e_corr: float = 0.4, 
                        rw_eps: float = 1e-36, 
                        rw_scalelim: float = 1.5,
                        grw_min_infeasible: int = 20, 
                        grw_interval: int = 20, 
                        grw_mod: float = 2.0):
        """
        Set parameters specific to the SDMM algorithm.

        Parameters
        ----------
        rw_interval : int, optional
            Interval for reweighting operations. 
            Defaults to 8.
        rw_e_corr : float, optional
            Error correction parameter for reweighting. 
            Defaults to 0.4.
        rw_eps : float, optional
            Small epsilon value for numerical stability in reweighting. 
            Defaults to 1e-36.
        rw_scalelim : float, optional
            Scale limit for reweighting operations. 
            Defaults to 1.5.
        grw_min_infeasible : int, optional
            Minimum number of infeasible iterations before adaptive reweighting is allowed. 
            Defaults to 20.
        grw_interval : int, optional
            Interval for adaptive reweighting checks. 
            Defaults to 20.
        grw_mod : float, optional
            Modification factor for adaptive reweighting. 
            Defaults to 2.0.
        """
        self.c_solver.set_sdmm_params(rw_interval, rw_e_corr, rw_eps, rw_scalelim,
                                     grw_min_infeasible, grw_interval, grw_mod)

    def solve(self, GroptParams gparams):
        """
        Run the SDMM solver.

        Uses the problem definitions set during the initialization of the solver.
        """
        self.c_solver.solve(gparams.c_gparams)

# ---------------------------------------------------------
# gropt_utils
# ---------------------------------------------------------
def set_verbose(level:int = 4, 
                mode: str | None = None):
    """
    Set the verbosity level for the C++ gropt library.

    The level to mode mapping is as follows:
    0: Off
    1: Critical
    2: Error
    3: Warning
    4: Info
    5: Debug
    6: Trace

    Parameters
    ----------  
    level : int
        The verbosity level to set.
    mode : str, optional
        The mode to set the verbosity level to. 
        Overrides the level if provided.
    """

    if mode is not None:
        if mode == 'off':
            level = 0
        elif mode == 'critical':
            level = 1
        elif mode == 'error':
            level = 2
        elif mode == 'warning':
            level = 3
        elif mode == 'info':
            level = 4
        elif mode == 'debug':
            level = 5
        elif mode == 'trace':
            level = 6
        else:
            raise ValueError(f"Unknown verbosity mode: {mode}")
    
    c_gropt.set_verbose(level)

def get_SAFE(G: np.ndarray, 
             dt: float, 
             true_safe: bool = True, 
             new_first_axis: int = 0, 
             demo_params: bool = True, 
             safe_params: dict = None) -> np.ndarray:

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

def get_SAFE2(G: np.ndarray, 
             dt: float, 
             true_safe: bool = True, 
             new_first_axis: int = 0, 
             demo_params: bool = True, 
             safe_params: dict = None) -> np.ndarray:

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

        c_gropt.get_SAFE2(N, Naxis, dt, &G_view[0], true_safe, new_first_axis, False,
                       &tau1_view[0], &tau2_view[0], &tau3_view[0],
                       &a1_view[0], &a2_view[0], &a3_view[0],
                       &stim_limit_view[0], &g_scale_view[0],
                       &out, out_size)
    else:
        c_gropt.get_SAFE2(N, Naxis, dt, &G_view[0], true_safe, new_first_axis, demo_params,
                        _unused, _unused, _unused, 
                        _unused, _unused, _unused,
                        _unused, _unused,
                        &out, out_size)

    
    return np.asarray(<cnp.float64_t[:out_size]> out)