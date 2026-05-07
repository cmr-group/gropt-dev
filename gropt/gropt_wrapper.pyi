from collections.abc import Callable
import enum
from typing import Annotated, overload

import numpy
from numpy.typing import NDArray


__build_date__: str = 'May  6 2026 20:08:25'

def set_log_level(level: int) -> None:
    """
    Set the log level for the C++ gropt library.

    Uses spdlog/Python logging convention: lower = more verbose.

    Level mapping: 0=Trace, 1=Debug, 2=Info, 3=Warning, 4=Error, 5=Critical, 6=Off.

    Parameters
    ----------
    level : int
        Log level (0–6).
    """

def set_log_callback(arg: Callable, /) -> None: ...

class SolveResult:
    """
    Result from a GrOpt solve operation.

    Attributes
    ----------
    X : np.ndarray
        The optimized gradient waveform.
    converged : bool
        Whether all constraints were satisfied.
    n_iter : int
        Number of outer SDMM iterations.
    n_feval : int
        Total number of inner linear solver iterations.
    """

    def __init__(self) -> None: ...

    @property
    def X(self) -> Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')]: ...

    @X.setter
    def X(self, arg: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')], /) -> None: ...

    @property
    def converged(self) -> bool: ...

    @converged.setter
    def converged(self, arg: bool, /) -> None: ...

    @property
    def n_iter(self) -> int: ...

    @n_iter.setter
    def n_iter(self, arg: int, /) -> None: ...

    @property
    def n_feval(self) -> int: ...

    @n_feval.setter
    def n_feval(self, arg: int, /) -> None: ...

    @property
    def dt(self) -> float: ...

    @dt.setter
    def dt(self, arg: float, /) -> None: ...

    @property
    def bvalue(self) -> float: ...

    @bvalue.setter
    def bvalue(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

class GroptParams:
    """
    Problem definition for GrOpt gradient optimization.

    This class holds the waveform dimensions, constraints, and objectives
    that define a gradient optimization problem.
    """

    def __init__(self) -> None: ...

    @property
    def N(self) -> int:
        """Number of time points per axis."""

    @N.setter
    def N(self, arg: int, /) -> None: ...

    @property
    def Naxis(self) -> int:
        """Number of gradient axes."""

    @Naxis.setter
    def Naxis(self, arg: int, /) -> None: ...

    @property
    def dt(self) -> float:
        """Raster time in seconds."""

    @dt.setter
    def dt(self, arg: float, /) -> None: ...

    def vec_init_simple(self, N: int = -1, Naxis: int = -1, first_val: float = 0.0, last_val: float = 0.0) -> None:
        """
        Initialize the set_vec and inv_vec settings.

        Parameters
        ----------
        N : int, optional
            Number of points in a single axis. Negative values use existing value.
        Naxis : int, optional
            Number of axes. Negative values use existing value.
        first_val : float, optional
            Fixed value for the first point [mT/m].
        last_val : float, optional
            Fixed value for the last point [mT/m].
        """

    def diff_init(self, dt: float = 0.0004, TE: float = 0.08, T_90: float = 0.003, T_180: float = 0.005, T_readout: float = 0.016) -> None:
        """
        Initialize diffusion sequence parameters.

        Parameters
        ----------
        dt : float, optional
            Raster time in seconds.
        TE : float, optional
            Echo time in seconds.
        T_90 : float, optional
            Duration of excitation RF pulse in seconds (time from excitation,
            so half of full RF pulse duration).
        T_180 : float, optional
            Duration of refocusing RF pulse in seconds.
        T_readout : float, optional
            Time to TE of the readout in seconds.
        """

    def diff_init_deadtime(self, dt: float = 0.0004, TE: float = 0.08, T_90: float = 0.003, T_180: float = 0.005, T_readout: float = 0.016) -> None:
        """
        Initialize diffusion sequence parameters, but with forced deadtime to make "convnetional" waveforms.

        Parameters
        ----------
        dt : float, optional
            Raster time in seconds.
        TE : float, optional
            Echo time in seconds.
        T_90 : float, optional
            Duration of excitation RF pulse in seconds (time from excitation,
            so half of full RF pulse duration).
        T_180 : float, optional
            Duration of refocusing RF pulse in seconds.
        T_readout : float, optional
            Time to TE of the readout in seconds.
        """

    def diff_init_preencode(self, dt: float = 0.0004, TE: float = 0.08, T_90: float = 0.003, T_180: float = 0.005, T_readout: float = 0.016, T_pre: float = 0.0) -> int:
        """
        Initialize diffusion sequence parameters with a pre-encoding period.

        Parameters
        ----------
        dt : float, optional
            Raster time in seconds.
        TE : float, optional
            Echo time in seconds.
        T_90 : float, optional
            Duration of excitation RF pulse in seconds (time from excitation,
            so half of full RF pulse duration).
        T_180 : float, optional
            Duration of refocusing RF pulse in seconds.
        T_readout : float, optional
            Time to TE of the readout in seconds.
        T_pre : float, optional
            Duration of the pre-encoding period in seconds.

        Returns
        -------
        int
            Number of pre-encoding time points (N_pre).
        """

    @overload
    def setvec_X0(self, X0: Annotated[NDArray[numpy.float64], dict(shape=(None,))], set_others: bool = True) -> None:
        """
        Set the initial waveform guess (1D).

        Parameters
        ----------
        X0 : np.ndarray
            Initial guess for the waveform (warm start) [T/m].
        set_others : bool, optional
            If True, update inv_vec, set_vals, and fixer with standard values.
        """

    @overload
    def setvec_X0(self, X0: Annotated[NDArray[numpy.float64], dict(shape=(None, None))], set_others: bool = True) -> None:
        """
        Set the initial waveform guess (2D: Naxis x N).

        Parameters
        ----------
        X0 : np.ndarray
            Initial guess for the waveform (warm start) [T/m]. Shape (Naxis, N).
        set_others : bool, optional
            If True, update inv_vec, set_vals, and fixer with standard values.
        """

    @overload
    def setvec_set_vals(self, set_vals: Annotated[NDArray[numpy.float64], dict(shape=(None,))]) -> None:
        """
        Manually set the set_vals vector (1D), and update fixer from its NaN pattern.

        NaN entries mark a point as free, finite entries lock the waveform to that value
        at that index. The fixer mask is rebuilt automatically (1.0 = free, 0.0 = fixed),
        and X0 is updated at the fixed positions to match.

        Note: pdata.inv_vec is not touched — call vec_init_simple() (or one of the
        diff_init variants) first if you need it allocated.

        Parameters
        ----------
        set_vals : np.ndarray
            1D array of length N. NaN = free, finite value = fixed.
        """

    @overload
    def setvec_set_vals(self, set_vals: Annotated[NDArray[numpy.float64], dict(shape=(None, None))]) -> None:
        """
        Manually set the set_vals vector (2D: Naxis x N), and update fixer from its NaN pattern.

        NaN entries mark a point as free, finite entries lock the waveform to that value
        at that index. The fixer mask is rebuilt automatically, and X0 is updated at the
        fixed positions to match.

        Parameters
        ----------
        set_vals : np.ndarray
            2D array of shape (Naxis, N). NaN = free, finite value = fixed.
        """

    def set_ils_solver(self, ils_method: str = 'CG') -> None:
        """
        Set the indirect solver method.

        Parameters
        ----------
        ils_method : str
            Solver method: 'CG', 'NLCG', or 'BiCGstabl' (case-sensitive).
        """

    def add_gmax(self, gmax: float = 0.03, rot_variant: bool = True, weight_mod: float = 1.0) -> None:
        """
        Add a maximum gradient amplitude constraint.

        Parameters
        ----------
        gmax : float, optional
            Maximum allowed gradient magnitude [T/m].
        rot_variant : bool, optional
            If True, use rotationally invariant formulation.
        weight_mod : float, optional
            Weighting factor for this constraint.
        """

    def add_concomitant(self, start_idx: int = 0, rot_variant: bool = True, weight_mod: float = 1.0) -> None:
        """
        Add a concomitant constraint.

        Parameters
        ----------
        start_idx : int, optional
            Starting index for the constraint (-1 = beginning).
        rot_variant : bool, optional
            If True, use rotationally invariant formulation.
        weight_mod : float, optional
            Weighting factor for this constraint.
        """

    def add_smax(self, smax: float = 80.0, rot_variant: bool = True, weight_mod: float = 1.0) -> None:
        """
        Add a maximum gradient slew rate constraint.

        Parameters
        ----------
        smax : float, optional
            Maximum allowed gradient slew rate [T/m/s].
        rot_variant : bool, optional
            If True, use rotationally invariant formulation.
        weight_mod : float, optional
            Weighting factor for this constraint.
        """

    def add_moment(self, order: float = 0, target: float = 0.0, tol: float = 1e-06, units: str = 'mT*ms/m', axis: int = 0, start_idx: int = -1, stop_idx: int = -1, ref_idx: int = 0, weight_mod: float = 1.0) -> None:
        """
        Add a moment constraint.

        Parameters
        ----------
        order : int, optional
            Moment order.
        target : float, optional
            Target moment value.
        tol : float, optional
            Tolerance for satisfying the constraint.
        units : str, optional
            Units: 'mT*ms/m', 'T*s/m', 'rad*s/m', or 's/m'.
        axis : int, optional
            Axis for moment calculation.
        start_idx : int, optional
            Starting index (-1 = beginning).
        stop_idx : int, optional
            Stopping index (-1 = end).
        ref_idx : int, optional
            Reference index (t=0 for moment calculations).
        weight_mod : float, optional
            Weighting factor for this constraint.
        """

    def add_SAFE(self, stim_thresh: float = 1.0, new_first_axis: int = 0, demo_params: bool = True, safe_params: object | None = None, weight_mod: float = 1.0) -> None:
        """
        Add a SAFE (PNS) constraint.

        Parameters
        ----------
        stim_thresh : float, optional
            Stimulus threshold for the SAFE constraint.
        new_first_axis : int, optional
            Swap the first axis of the SAFE parameters.
        demo_params : bool, optional
            Whether to use demo parameters. Ignored if safe_params is provided.
        safe_params : dict, optional
            Dictionary of SAFE parameters (see gropt.readasc).
        weight_mod : float, optional
            Weighting factor for this constraint.
        """

    def add_SAFE_vec(self, stim_thresh_vec: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')], new_first_axis: int = 0, demo_params: bool = True, safe_params: object | None = None, weight_mod: float = 1.0) -> None:
        """
        Add a SAFE constraint with a vector stimulation limit.

        Parameters
        ----------
        stim_thresh_vec : np.ndarray
            Vector of stimulus thresholds.
        new_first_axis : int, optional
            Swap the first axis of the SAFE parameters.
        demo_params : bool, optional
            Whether to use demo parameters.
        safe_params : dict, optional
            Dictionary of SAFE parameters (see gropt.readasc).
        weight_mod : float, optional
            Weighting factor for this constraint.
        """

    def add_eddy(self, lam: object, tol: float = 0.0001, weight_mod: float = 1.0) -> None:
        """
        Add an eddy current constraint.

        Parameters
        ----------
        lam : float or np.ndarray
            Time constant(s) for eddy currents [seconds]. A single float is treated
            as a one-element array.
        tol : float, optional
            Tolerance for the constraint.
        weight_mod : float, optional
            Weighting factor for this constraint.
        """

    def add_bvalue(self, target: float = 100.0, tol: float = 1.0, start_idx0: int = -1, stop_idx0: int = -1, weight_mod: float = 1.0, mode: object = 2, max_scale: float = 1.01) -> None:
        """
        Add a b-value constraint.

        Parameters
        ----------
        target : float, optional
            Target b-value.
        tol : float, optional
            Tolerance for the b-value constraint.
        start_idx0 : int, optional
            Starting index (-1 = full waveform).
        stop_idx0 : int, optional
            Stopping index (-1 = full waveform).
        weight_mod : float, optional
            Weighting factor.
        mode : int or str, optional
            'setval'/1 = set b-value, 'minval'/2 = minimum b-value (default),
            'minval_max'/3 = minimum b-value with scaling.
        max_scale : float, optional
            Scale factor when mode=3.
        """

    def add_TV(self, tv_lam: float = 0.0, weight_mod: float = 1.0) -> None:
        """
        Add total variation regularization.

        Parameters
        ----------
        tv_lam : float, optional
            Regularization strength (must be > 0 to have effect).
        weight_mod : float, optional
            Weighting factor.
        """

    def add_obj_identity(self, weight_mod: float = 1.0) -> None:
        """
        Add an identity (L2 norm) objective.

        Parameters
        ----------
        weight_mod : float, optional
            Weighting factor.
        """

    def prepare(self) -> None:
        """
        Prepare the problem for solving.

        Allocates vectors and sets up initial optimization variables.
        Automatically called by solve() if not already done.
        """

    def print_op_details(self) -> None:
        """
        Print details of all operators.

        This function prints information about each operator in the problem, including their types and parameters.
        """

    def reset_op_weights(self) -> None:
        """
        Reset all operator weights and spectral norms to 1.0.

        Sets weight_mod, spec_norm, and spec_norm2 to 1.0 on every
        operator in all_op and all_obj.
        """

class Solver:
    @property
    def extra_debug(self) -> bool:
        """
        If True, populate debug_solver with per-iteration history after solve().
        """

    @extra_debug.setter
    def extra_debug(self, arg: bool, /) -> None: ...

    def get_debug(self) -> dict:
        """
        Return debug history as a dict of lists of numpy arrays.

        Only populated when extra_debug=True before solve().

        Returns
        -------
        dict with keys: hist_X, hist_Ax, hist_z, hist_y, hist_Aty
            Each value is a list of 1-D numpy arrays, one per logged iteration.
        """

    def set_general_params(self, min_iter: int = 1, max_iter: int = 2000, log_interval: int = 20, gamma_x: float = 1.6, max_feval: int = 12000, extra_iters: int = 0) -> None:
        """
        Set general solver parameters.

        Parameters
        ----------
        min_iter : int, optional
            Minimum iterations.
        max_iter : int, optional
            Maximum iterations.
        log_interval : int, optional
            Logging interval (only visible with verbose logging).
        gamma_x : float, optional
            Relaxation parameter for updates.
        max_feval : int, optional
            Maximum total function evaluations.
        extra_iters : int, optional
            Number of extra iterations to run after solution is found.
        """

    def set_ils_params(self, ils_tol: float = 0.001, ils_max_iter: int = 20, ils_min_iter: int = 2, ils_sigma: float = 0.0001, ils_tik_lam: float = 0.0) -> None:
        """
        Set indirect linear solver parameters.

        Parameters
        ----------
        ils_tol : float, optional
            Relative tolerance for the inner solver.
        ils_max_iter : int, optional
            Maximum ILS iterations per outer iteration.
        ils_min_iter : int, optional
            Minimum ILS iterations.
        ils_sigma : float, optional
            ADMM penalty parameter.
        ils_tik_lam : float, optional
            Tikhonov regularization parameter.
        """

class SolverGroptSDMM(Solver):
    """SDMM solver for GrOpt gradient optimization problems."""

    def __init__(self) -> None: ...

    def set_sdmm_params(self, rw_interval: int = 8, rw_e_corr: float = 0.4, rw_eps: float = 1e-36, rw_scalelim: float = 1.5, grw_min_infeasible: int = 20, grw_interval: int = 20, grw_mod: float = 2.0) -> None:
        """
        Set SDMM-specific parameters.

        Parameters
        ----------
        rw_interval : int, optional
            Interval for reweighting operations.
        rw_e_corr : float, optional
            Error correction for reweighting.
        rw_eps : float, optional
            Epsilon for numerical stability in reweighting.
        rw_scalelim : float, optional
            Scale limit for reweighting.
        grw_min_infeasible : int, optional
            Minimum infeasible iterations before adaptive reweighting.
        grw_interval : int, optional
            Interval for adaptive reweighting checks.
        grw_mod : float, optional
            Modification factor for adaptive reweighting.
        """

    def solve(self, gparams: GroptParams) -> SolveResult:
        """
        Run the SDMM solver.

        Parameters
        ----------
        gparams : GroptParams
            The problem definition.

        Returns
        -------
        SolveResult
            The optimization result containing the waveform and convergence info.
        """

class SolverOSQP(Solver):
    """OSQP solver for GrOpt gradient optimization problems."""

    def __init__(self) -> None: ...

    def solve(self, gparams: GroptParams) -> SolveResult:
        """
        Run the OSQP solver.

        Parameters
        ----------
        gparams : GroptParams
            The problem definition.

        Returns
        -------
        SolveResult
            The optimization result containing the waveform and convergence info.
        """

def solve(params: GroptParams, min_iter: int = 1, max_iter: int = 2000, log_interval: int = 20, gamma_x: float = 1.6, max_feval: int = 12000, extra_iters: int = 0, ils_tol: float = 0.001, ils_max_iter: int = 20, ils_min_iter: int = 2, ils_sigma: float = 0.0001, ils_tik_lam: float = 0.0) -> SolveResult:
    """
    Convenience function to solve a GrOpt problem.

    Parameters
    ----------
    params : GroptParams
        The problem definition.
    min_iter, max_iter, log_interval, gamma_x, max_feval
        General solver parameters (see SolverGroptSDMM.set_general_params).
    ils_tol, ils_max_iter, ils_min_iter, ils_sigma, ils_tik_lam
        ILS parameters (see SolverGroptSDMM.set_ils_params).

    Returns
    -------
    SolveResult
        The optimization result.
    """

class NormType(enum.Enum):
    L2 = 0

    Inf = 1

def estimate_row_col_norms(gparams: GroptParams, n_reps: int = 10, norm_type: NormType = NormType.Inf) -> tuple[Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')], Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')]]:
    """
    Estimate row and column norms of the operator matrix.

    Parameters
    ----------
    gparams : GroptParams
        The problem definition (must have operators added).
    n_reps : int, optional
        Number of random vector repetitions for the estimate.
    norm_type : NormType, optional
        NormType.Inf (default) or NormType.L2.

    Returns
    -------
    row_norms : np.ndarray
        Estimated norm for each constraint row.
    col_norms : np.ndarray
        Estimated norm for each variable column.
    """

def get_eq_vecs(gparams: GroptParams) -> tuple[Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')], Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')]]:
    """
    Get the current equilibration vectors from all operators.

    Returns the accumulated eq_rows and eq_cols stored on each operator
    after a call to equilibrate().

    Parameters
    ----------
    gparams : GroptParams
        The problem definition (must have operators prepared and equilibrated).

    Returns
    -------
    row_norms : np.ndarray
        Accumulated row equilibration vector.
    col_norms : np.ndarray
        Accumulated column equilibration vector.
    """

def equilibrate(gparams: GroptParams, n_iter: int = 5, n_reps: int = 10) -> None:
    """
    Ruiz equilibration of the operator matrix.

    Iteratively scales operator rows and columns so that the maximum
    absolute value in each row and column approaches 1, improving
    solver conditioning.

    Parameters
    ----------
    gparams : GroptParams
        The problem definition (operators will be modified in place).
    n_iter : int, optional
        Number of equilibration iterations.
    n_reps : int, optional
        Number of random vector repetitions per norm estimate.
    """

def rescale_eq_vecs(gparams: GroptParams, row_scale: float, col_scale: float) -> None:
    """
    Rescale equilibration vectors by scalar factors.

    Multiplies all operator eq_rows by row_scale and eq_cols by col_scale.

    Parameters
    ----------
    gparams : GroptParams
        The problem definition.
    row_scale : float
        Scale factor applied to all row equilibration vectors.
    col_scale : float
        Scale factor applied to all column equilibration vectors.
    """

def estimate_spec_norm(gparams: GroptParams, n_iters: int = 20) -> float:
    """
    Estimate the spectral norm of the operator matrix via power iteration.

    Parameters
    ----------
    gparams : GroptParams
        The problem definition (must have operators prepared).
    n_iters : int, optional
        Number of power iterations.

    Returns
    -------
    float
        Estimated spectral norm.
    """

def estimate_individual_spec_norm(gparams: GroptParams, n_iters: int = 20, op_idx: int = 0) -> float:
    """
    Estimate the spectral norm of a single operator matrix via power iteration.

    Parameters
    ----------
    gparams : GroptParams
        The problem definition (must have operators prepared).
    n_iters : int, optional
        Number of power iterations.
    op_idx : int, optional
        Index of the operator for which to estimate the spectral norm.

    Returns
    -------
    float
        Estimated spectral norm.
    """

def get_SAFE(G: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')], dt: float, true_safe: bool = True, new_first_axis: int = 0, demo_params: bool = True, safe_params: object | None = None) -> Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')]:
    """
    Compute the SAFE (PNS) response for a gradient waveform.

    Parameters
    ----------
    G : np.ndarray
        Gradient waveform.
    dt : float
        Raster time in seconds.
    true_safe : bool, optional
        Use true SAFE model.
    new_first_axis : int, optional
        Swap the first axis of SAFE parameters.
    demo_params : bool, optional
        Use demo parameters.
    safe_params : dict, optional
        Dictionary of SAFE parameters (see gropt.readasc).

    Returns
    -------
    np.ndarray
        SAFE response curve.
    """

def test_eigen_assertions(test_type: int) -> None:
    """
    Test that assertions are running in Eigen. 

    Each test is a different assertion that should be triggered in Eigen.  All tests are
    expected to pass in normal execution.

    Parameters
    ----------
    test_type : int
        The type of assertion test to run, (1, 2, or 3).
    """
