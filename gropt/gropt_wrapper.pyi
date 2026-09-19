"""GrOpt: Gradient Optimization for MRI"""

from collections.abc import Callable
import enum
from typing import Annotated, overload

import numpy
from numpy.typing import NDArray


__build_date__: str = 'Sep 18 2026 21:36:42'

def set_log_level(level: int) -> None:
    """
    Set the log level for the C++ gropt library (lower = more verbose).

    Parameters
    ----------
    level : int
        0=Trace, 1=Debug, 2=Info, 3=Warning, 4=Error, 5=Critical, 6=Off.
    """

def set_log_callback(arg: Callable, /) -> None:
    """
    Route C++ log messages to a Python callable.

    Replaces all sinks of the C++ default logger with one that calls fn(level, message),
    where level is an int as in set_log_level and message is a str. Prefer
    gropt.setup_logging(), which also releases the callback at exit.
    """

def clear_log_callback() -> None:
    """
    Remove all C++ log sinks, releasing the Python callback; setup_logging registers this at exit.
    """

class SolveResult:
    """
    Result from a GrOpt solve operation.

    Attributes
    ----------
    X : np.ndarray
        The optimized gradient waveform [T/m], flat, length Naxis*N (axis-major).
    converged : bool
        True if every constraint is feasible for the returned waveform.
    n_iter : int
        Number of outer iterations.
    n_feval : int
        Total number of inner linear-solver iterations.
    dt : float
        Raster time of the waveform [s].
    bvalue : float
        b-value of the returned waveform [s/mm^2] if the problem has a b-value term
        (constraint or objective), else 0.
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
    """Problem definition: waveform layout, constraints, and objectives."""

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
        """Raster time [s]."""

    @dt.setter
    def dt(self, arg: float, /) -> None: ...

    @property
    def normalize_obj(self) -> bool:
        """
        If True, scale each linearized objective's pull (e.g. add_bvalue with as_objective=True) to unit norm times its weight, so weight_mod acts as a step size. Default False.
        """

    @normalize_obj.setter
    def normalize_obj(self, arg: bool, /) -> None: ...

    @property
    def eq_proj_solver(self) -> EqProjSolver:
        """
        Factorization for the project=True equality projection. EqProjSolver.LDLT (default): fast, but fails on singular or collinear rows. EqProjSolver.COD: rank-revealing, handles near-collinear rows (e.g. several eddy time constants).
        """

    @eq_proj_solver.setter
    def eq_proj_solver(self, arg: EqProjSolver, /) -> None: ...

    @property
    def eq_proj_rcond(self) -> float:
        """
        Rank tolerance for EqProjSolver.COD: singular values below eq_proj_rcond times the largest are dropped; <= 0 uses Eigen's default. Ignored by LDLT. Default 1e-10.
        """

    @eq_proj_rcond.setter
    def eq_proj_rcond(self, arg: float, /) -> None: ...

    @property
    def safe_eps(self) -> float:
        """
        Softabs smoothing [T/m/s] for SAFE |.|: sqrt(v^2 + eps^2), slightly conservative; 0 = exact (default). Set before add_SAFE/add_SAFE_vec; try ~1% of smax (e.g. 1-5 for smax = 200).
        """

    @safe_eps.setter
    def safe_eps(self, arg: float, /) -> None: ...

    def vec_init_simple(self, N: int = -1, Naxis: int = -1, first_val: float = 0.0, last_val: float = 0.0) -> None:
        """
        Initialize a simple (non-diffusion) problem layout.

        Sets inv_vec to +1, fixes the first and last points of every axis to first_val and
        last_val (all other points free), and sets X0 to 0.01 at the free points.

        Parameters
        ----------
        N : int, optional
            Number of points per axis; <= 0 keeps the current value.
        Naxis : int, optional
            Number of axes; <= 0 keeps the current value.
        first_val : float, optional
            Fixed value for the first point [T/m].
        last_val : float, optional
            Fixed value for the last point [T/m].
        """

    def diff_init(self, dt: float = 0.0004, TE: float = 0.08, T_90: float = 0.003, T_180: float = 0.005, T_readout: float = 0.016) -> None:
        """
        Initialize a single-axis spin-echo diffusion layout (N, inv_vec, fixed RF blocks, X0).

        Parameters
        ----------
        dt : float, optional
            Raster time [s].
        TE : float, optional
            Echo time [s].
        T_90 : float, optional
            Excitation RF time after t = 0 [s] (half the full pulse duration).
        T_180 : float, optional
            Refocusing RF pulse duration [s], centered at TE/2.
        T_readout : float, optional
            Readout time before TE [s]; the waveform ends at TE - T_readout.
        """

    def diff_init_deadtime(self, dt: float = 0.0004, TE: float = 0.08, T_90: float = 0.003, T_180: float = 0.005, T_readout: float = 0.016) -> None:
        """
        Like diff_init, with dead time for "conventional" waveforms.

        The part of the pre-180 period that exceeds the post-180 period is fixed to
        zero, so both sides have the same free duration.

        Parameters
        ----------
        dt : float, optional
            Raster time [s].
        TE : float, optional
            Echo time [s].
        T_90 : float, optional
            Excitation RF time after t = 0 [s] (half the full pulse duration).
        T_180 : float, optional
            Refocusing RF pulse duration [s], centered at TE/2.
        T_readout : float, optional
            Readout time before TE [s]; the waveform ends at TE - T_readout.
        """

    def diff_init_preencode(self, dt: float = 0.0004, TE: float = 0.08, T_90: float = 0.003, T_180: float = 0.005, T_readout: float = 0.016, T_pre: float = 0.0) -> int:
        """
        Like diff_init, with a free pre-encoding period prepended.

        Parameters
        ----------
        dt : float, optional
            Raster time [s].
        TE : float, optional
            Echo time [s].
        T_90 : float, optional
            Excitation RF time after t = 0 [s] (half the full pulse duration).
        T_180 : float, optional
            Refocusing RF pulse duration [s], centered at TE/2.
        T_readout : float, optional
            Readout time before TE [s]; the waveform ends at TE - T_readout.
        T_pre : float, optional
            Duration of the pre-encoding period [s].

        Returns
        -------
        int
            Number of pre-encoding time points (N_pre).
        """

    @overload
    def setvec_X0(self, X0: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C', writable=False)], set_others: bool = True) -> None:
        """
        Set the initial waveform guess (1D: single axis).

        Parameters
        ----------
        X0 : np.ndarray
            Initial guess [T/m], length N. Sets Naxis = 1.
        set_others : bool, optional
            If True, reset inv_vec to +1 and fix each axis's first and last points to
            their X0 values (all others free). Use False to keep a layout from diff_init.
        """

    @overload
    def setvec_X0(self, X0: Annotated[NDArray[numpy.float64], dict(shape=(None, None), order='C', writable=False)], set_others: bool = True) -> None:
        """
        Set the initial waveform guess (2D: Naxis x N).

        Parameters
        ----------
        X0 : np.ndarray
            Initial guess [T/m], shape (Naxis, N).
        set_others : bool, optional
            If True, reset inv_vec to +1 and fix each axis's first and last points to
            their X0 values (all others free). Use False to keep a layout from diff_init.
        """

    @overload
    def setvec_set_vals(self, set_vals: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C', writable=False)]) -> None:
        """
        Set the fixed-value vector set_vals (1D: single axis) and rebuild the fixer mask.

        X0 is set to the fixed values at fixed points. inv_vec is not changed; call
        vec_init_simple() or a diff_init variant first if it is not allocated.

        Parameters
        ----------
        set_vals : np.ndarray
            Length N; sets Naxis = 1. NaN = free, finite value = fixed to that value.
        """

    @overload
    def setvec_set_vals(self, set_vals: Annotated[NDArray[numpy.float64], dict(shape=(None, None), order='C', writable=False)]) -> None:
        """
        Set the fixed-value vector set_vals (2D: Naxis x N) and rebuild the fixer mask.

        X0 is set to the fixed values at fixed points. inv_vec is not changed.

        Parameters
        ----------
        set_vals : np.ndarray
            Shape (Naxis, N). NaN = free, finite value = fixed to that value.
        """

    def getvec_set_vals(self) -> Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')]:
        """
        Return a copy of set_vals (flat, length Naxis*N): NaN = free, finite = fixed.
        """

    def getvec_fixer(self) -> Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')]:
        """
        Return a copy of the fixer mask (flat, length Naxis*N): 1.0 = free, 0.0 = fixed.
        """

    def getvec_inv_vec(self) -> Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')]:
        """
        Return a copy of inv_vec (flat, length Naxis*N): the +/- inversion pattern (sign flip at the 180).
        """

    def getvec_X0(self) -> Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')]:
        """Return a copy of the current initial-guess X0 (flat, length Naxis*N)."""

    def set_ils_solver(self, ils_method: str = 'CG') -> None:
        """
        Set the inner (indirect) linear solver.

        Only CG applies the equality projection used by project=True constraints; with
        NLCG or BiCGstabl they are enforced only by reproject_iterate (slower).

        Parameters
        ----------
        ils_method : str
            Solver method: 'CG' (default), 'NLCG', or 'BiCGstabl' (case-sensitive).
        """

    def add_gmax(self, gmax: float = 0.03, rot_variant: bool = True, weight_mod: float = 1.0) -> None:
        """
        Add a maximum gradient amplitude constraint.

        Parameters
        ----------
        gmax : float, optional
            Maximum allowed gradient magnitude [T/m].
        rot_variant : bool, optional
            If True (default), limit each axis independently. If False, limit the
            gradient magnitude across axes, which is rotationally invariant.
        weight_mod : float, optional
            Weighting factor for this constraint.
        """

    def add_gmax_vec(self, gmax_vec: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')], rot_variant: bool = True, weight_mod: float = 1.0) -> None:
        """
        Add a per-location gradient amplitude constraint.

        Parameters
        ----------
        gmax_vec : np.ndarray
            Gradient amplitude limit at each location [T/m]. Length N is broadcast
            across all axes; length Naxis*N sets a per-axis limit.
        rot_variant : bool, optional
            Must be True; a rotationally invariant vector limit is not supported.
        weight_mod : float, optional
            Weighting factor for this constraint.
        """

    def add_concomitant(self, start_idx: int = 0, rot_variant: bool = True, weight_mod: float = 1.0, tol0: float = 0.1, target: float = 1.0, project: bool = False) -> None:
        """
        Add a concomitant (pre/post-180 energy balance) constraint.

        Constrains the ratio pos/neg of the gradient energy sum(g^2 * dt) before (pos)
        and after (neg) the 180 to target +/- tol0. The constraint is nonconvex: a soft
        ADMM constraint by default, or a relinearized equality projection with
        project=True. The prox is the exact projection onto the ratio band.

        Parameters
        ----------
        start_idx : int, optional
            Index where the pre/post energy sums start (0 = beginning).
        rot_variant : bool, optional
            Ignored; the balance is always computed across all axes.
        weight_mod : float, optional
            Weighting factor for this constraint.
        tol0 : float, optional
            Tolerance on the ratio: feasible when |pos/neg - target| <= tol0.
        target : float, optional
            Target ratio pos/neg (default 1.0 = balanced). With project=True the
            linearized pos - target*neg = 0 is enforced.
        project : bool, optional
            If True, enforce the linearized (SQP) balance via the equality projection
            each outer iteration; weight_mod is unused and tol0 only sets the
            feasibility check.
        """

    def add_smax(self, smax: float = 80.0, rot_variant: bool = True, weight_mod: float = 1.0) -> None:
        """
        Add a maximum gradient slew rate constraint.

        Parameters
        ----------
        smax : float, optional
            Maximum allowed gradient slew rate [T/m/s].
        rot_variant : bool, optional
            If True (default), limit each axis independently. If False, limit the
            slew magnitude across axes, which is rotationally invariant.
        weight_mod : float, optional
            Weighting factor for this constraint.
        """

    def add_smax_vec(self, smax_vec: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')], rot_variant: bool = True, weight_mod: float = 1.0) -> None:
        """
        Add a per-location slew rate constraint.

        Parameters
        ----------
        smax_vec : np.ndarray
            Slew-rate limit at each location [T/m/s]. Length (N-1) is broadcast
            across all axes; length Naxis*(N-1) sets a per-axis limit.
        rot_variant : bool, optional
            Must be True; a rotationally invariant vector limit is not supported.
        weight_mod : float, optional
            Weighting factor for this constraint.
        """

    def add_moment(self, order: float = 0, target: float = 0.0, tol: float = 1e-06, units: str = 'mT*ms/m', axis: int = 0, start_idx: int = -1, stop_idx: int = -1, ref_idx: int = 0, weight_mod: float = 1.0, project: bool = False, absolute_tol: bool = False) -> None:
        """
        Add a moment constraint.

        Parameters
        ----------
        order : int, optional
            Moment order (0, 1, 2, ...).
        target : float, optional
            Target moment value, in `units`.
        tol : float, optional
            Tolerance in `units` at order 0. For order k it is scaled by ||A_k|| / ||A_0||
            (moment-row norms over the window) unless absolute_tol=True. With
            project=True it only sets the feasibility check.
        units : str, optional
            'mT*ms/m', 'T*s/m', 'rad*s/m', or 's/m'.
        axis : int, optional
            Axis index.
        start_idx : int, optional
            Starting index (-1 = beginning).
        stop_idx : int, optional
            Stopping index (-1 = end).
        ref_idx : int, optional
            Index taken as t = 0 for the moment.
        weight_mod : float, optional
            Weighting factor for this constraint.
        project : bool, optional
            If True, enforce the moment exactly by null-space projection instead of an
            ADMM penalty.
        absolute_tol : bool, optional
            If True, tol is in this order's own units instead of scaled from M0 (e.g.
            for a nonzero M2 target). Default False.
        """

    def add_SAFE(self, stim_thresh: float = 1.0, new_first_axis: int = 0, demo_params: bool = True, safe_params: object | None = None, weight_mod: float = 1.0) -> None:
        """
        Add a SAFE (PNS) constraint.

        Parameters
        ----------
        stim_thresh : float, optional
            Stimulation limit as a fraction of the SAFE threshold (1.0 = 100%).
        new_first_axis : int, optional
            Use the SAFE parameters of this axis (0, 1, 2) for the first gradient axis
            (swapped with axis 0).
        demo_params : bool, optional
            Use the built-in demo SAFE parameters. Must be True if safe_params is None;
            ignored if safe_params is provided.
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
            Per-sample stimulation limits as fractions of the SAFE threshold, length
            Naxis*N. Any other length falls back to 1.0 everywhere, with a warning.
        new_first_axis : int, optional
            Use the SAFE parameters of this axis (0, 1, 2) for the first gradient axis
            (swapped with axis 0).
        demo_params : bool, optional
            Use the built-in demo SAFE parameters. Must be True if safe_params is None;
            ignored if safe_params is provided.
        safe_params : dict, optional
            Dictionary of SAFE parameters (see gropt.readasc).
        weight_mod : float, optional
            Weighting factor for this constraint.
        """

    def add_eddy(self, lam: object, tol: float = 0.0001, weight_mod: float = 1.0, project: bool = False) -> None:
        """
        Add an eddy-current constraint on the residual eddy term at the end of the waveform.

        The target is 0: a soft box |eddy| <= tol per time constant by default, or held
        at 0 by the equality projection (like add_moment) with project=True.

        Parameters
        ----------
        lam : float or np.ndarray
            Time constant(s) [s]. One row is added per (axis, time constant).
        tol : float, optional
            Box half-width on the residual eddy term. With project=True it only sets the
            feasibility check.
        weight_mod : float, optional
            Weighting factor for this constraint (project=False only).
        project : bool, optional
            If True, enforce eddy == 0 by null-space projection. Keep the number of time
            constants small, since many projected equalities (plus moments/concomitant)
            can over-constrain the free points. Close time constants give nearly
            collinear rows; use eq_proj_solver = EqProjSolver.COD in that case.
        """

    def add_bvalue(self, target: float = 100.0, tol: float = 1.0, start_idx0: int = -1, stop_idx0: int = -1, weight_mod: float = 1.0, mode: object = 2, max_scale: float = 1.01, as_objective: bool = False, linearize: bool = True) -> None:
        """
        Add a b-value term (constraint by default, or a maximization objective).

        Parameters
        ----------
        target : float, optional
            Target b-value [s/mm^2]. Constraint only.
        tol : float, optional
            Tolerance on the b-value [s/mm^2], used by mode='setval'.
        start_idx0 : int, optional
            Starting index (-1 = full waveform).
        stop_idx0 : int, optional
            Stopping index (-1 = full waveform).
        weight_mod : float, optional
            Constraint: weighting factor. Objective: strength of the b-value pull
            (obj_weight = -weight_mod).
        mode : int or str, optional
            Constraint only. 'setval'/1: b = target +/- tol. 'minval'/2 (default):
            b >= target. 'minval_max'/3: b >= target, and each prox also scales the
            waveform up by max_scale to keep pushing b higher.
        max_scale : float, optional
            Per-iteration scale factor for mode='minval_max'.
        as_objective : bool, optional
            If True, maximize the b-value as an objective instead of constraining it.
        linearize : bool, optional
            Objective only. True (default, recommended): linearized into the RHS
            (DCA), keeping the CG system positive-definite. False (experimental):
            enters as curvature in the LHS, which can make the CG system indefinite.
        """

    def add_diff_basin(self, window_time: float, eps_factor: float, gmax: float, weight_mod: float = 1.0, same_sign: bool = False) -> None:
        """
        Add a diffusion basin-orientation constraint (single 180 in the middle).

        Requires the mean gradient in a short window on each side of the 180 to satisfy

            mean(g, pre window)  >= +eps
            mean(g, post window) <= -eps     (>= +eps with same_sign=True)

        with eps = eps_factor * gmax. The 180 is located from the inv_vec sign flip, and
        each window takes round(window_time / dt) free samples on its side, skipping the
        fixed RF block. Use it to pick which sign pattern b-value maximization converges
        to. Keep eps_factor small so the constraint is inactive at the optimum.

        Parameters
        ----------
        window_time : float
            Window width on each side of the 180 [s] (e.g. 1e-3).
        eps_factor : float
            Minimum |mean window gradient| as a fraction of gmax (e.g. 0.07).
        gmax : float
            Gradient limit [T/m], used to scale eps = eps_factor * gmax.
        weight_mod : float, optional
            Weighting factor for this constraint.
        same_sign : bool, optional
            False (default): opposite signs across the 180 (post <= -eps), the sign
            pattern of the global-max diffusion waveform. True: the same sign on both
            sides (post >= +eps).

        Notes
        -----
        Assumes a single axis and a single 180. The pre window is always >= +eps; keep
        any X0 seed orientation consistent with that.
        """

    def add_TV(self, tv_lam: float = 0.0, weight_mod: float = 1.0, order: int = 1) -> None:
        """
        Add a total-variation (L1) penalty on a finite difference of the gradient.

        Minimizes tv_lam * ||D g||_1, where D g is the slew [T/m/s] (order 1) or jerk
        [T/m/s^2] (order 2). It is a penalty, not a constraint: it never affects
        feasibility or the choice of returned iterate, so with no other objective the
        result depends on min_iter. tv_lam is in physical units, so useful values are
        small and dt-dependent (e.g. 1e-9 to 1e-6 for order 2 at 400 us against a
        b-value objective; larger values saturate at the minimum-TV waveform).

        Parameters
        ----------
        tv_lam : float, optional
            Penalty weight; must be > 0 to have effect.
        weight_mod : float, optional
            Initial ADMM weight (rho). Affects convergence, not the penalty strength.
        order : int, optional
            1 = slew: sparse slew, blocky gradients; fights bang-bang ramps.
            2 = jerk: piecewise-constant slew; removes slew jitter without fighting
            the ramps.
        """

    def add_obj_identity(self, weight_mod: float = 1.0) -> None:
        """
        Add an identity (L2 norm) objective that minimizes ||g||^2.

        Parameters
        ----------
        weight_mod : float, optional
            Weighting factor.
        """

    def prepare(self) -> None:
        """
        Initialize all operators for the current layout.

        solve() calls this automatically if N changed or operators were added since
        the last prepare().
        """

    def print_op_details(self) -> None:
        """
        Log each constraint and objective operator's type and parameters (at Info level).
        """

    def reset_op_weights(self) -> None:
        """
        Set spec_norm and spec_norm2 to 1.0 on every constraint and objective operator.

        weight_mod is not changed.
        """

    def get_op_names(self) -> list[str]:
        """
        Return the constraint-operator names, in all_op order.

        This is the per-operator order of the Solver.get_debug() histories, so the list
        can be used directly as plot labels.
        """

    def get_op_keys(self) -> list[str]:
        """
        Return each constraint operator's unique key ("<name>#<occurrence>"), in all_op order.

        These are the keys warm-start snapshots use to match operators between solves.
        """

class Solver:
    @property
    def extra_debug(self) -> bool:
        """
        If True, record per-iteration histories during solve(); read them with get_debug().
        """

    @extra_debug.setter
    def extra_debug(self, arg: bool, /) -> None: ...

    @property
    def min_iter(self) -> int:
        """Minimum outer iterations before feasible iterates are considered."""

    @min_iter.setter
    def min_iter(self, arg: int, /) -> None: ...

    @property
    def max_iter(self) -> int:
        """Maximum outer iterations."""

    @max_iter.setter
    def max_iter(self, arg: int, /) -> None: ...

    @property
    def log_interval(self) -> int:
        """Outer iterations between debug log prints."""

    @log_interval.setter
    def log_interval(self, arg: int, /) -> None: ...

    @property
    def gamma_x(self) -> float:
        """
        Over-relaxation factor for the outer primal update, X <- gamma_x*Xhat + (1-gamma_x)*X.
        """

    @gamma_x.setter
    def gamma_x(self, arg: float, /) -> None: ...

    @property
    def max_feval(self) -> int:
        """Maximum total inner (CG) iterations over the whole solve."""

    @max_feval.setter
    def max_feval(self, arg: int, /) -> None: ...

    @property
    def obj_patience(self) -> int:
        """
        Objective problems: stop after this many feasible iters with no objective improvement.
        """

    @obj_patience.setter
    def obj_patience(self, arg: int, /) -> None: ...

    @property
    def obj_rtol(self) -> float:
        """
        Relative objective-improvement threshold for the obj_patience plateau test.
        """

    @obj_rtol.setter
    def obj_rtol(self, arg: float, /) -> None: ...

    @property
    def ils_tol(self) -> float:
        """
        Inner-solver tolerance, relative to the warm-start residual (~0.1 for inexact solves, smaller for near-exact).
        """

    @ils_tol.setter
    def ils_tol(self, arg: float, /) -> None: ...

    @property
    def ils_max_iter(self) -> int:
        """Max inner-solver iterations per outer step."""

    @ils_max_iter.setter
    def ils_max_iter(self, arg: int, /) -> None: ...

    @property
    def ils_min_iter(self) -> int:
        """Min inner-solver iterations per outer step."""

    @ils_min_iter.setter
    def ils_min_iter(self, arg: int, /) -> None: ...

    @property
    def ils_sigma(self) -> float:
        """
        Proximal weight: adds sigma*I to the inner system, anchoring each inner solve to the current iterate.
        """

    @ils_sigma.setter
    def ils_sigma(self, arg: float, /) -> None: ...

    @property
    def ils_tik_lam(self) -> float:
        """
        Tikhonov weight: adds tik_lam*I to the inner-system matrix, shrinking toward zero.
        """

    @ils_tik_lam.setter
    def ils_tik_lam(self, arg: float, /) -> None: ...

    def get_debug(self) -> dict:
        """
        Return the debug histories of the last solve() as a dict.

        Only populated when extra_debug=True before solve(). Lists have one entry per
        outer iteration unless noted. Per-operator lists are ordered like
        GroptParams.get_op_names(), except hist_r_prim and hist_r_dual, which skip
        operators with project=True.

        Returns
        -------
        dict with keys:
            hist_X, hist_Ax, hist_z, hist_y, hist_Aty
                list of 1-D numpy arrays.
            hist_weight, hist_gamma, hist_r_feas, hist_feas
                list of per-operator lists: ADMM weight, relaxation gamma, relative
                infeasibility, and binary feasible flag.
            hist_r_prim, hist_r_dual
                list of per-operator lists: Boyd ADMM primal/dual residuals.
            hist_all_feas
                list of ints: 1 if all operators were feasible that iteration, else 0.
            hist_gamma_x
                list of floats: the outer relaxation factor gamma_x.
            hist_bvalue
                b-value of the iterate (empty list if there is no b-value term).
            best_feasible_iter
                int: iteration of the returned iterate, or -1 if none was feasible
                (e.g. hist_bvalue[best_feasible_iter] is the returned b-value).
            hist_cg_iter, hist_cg_rnorm0, hist_cg_rnorm, hist_cg_bnorm0
                inner-solver diagnostics, one entry per inner solve (trust-region
                re-solves add entries): iterations taken, initial/final residual norm,
                and ||b||. Only hist_cg_iter has a leading -1 placeholder, so
                hist_cg_iter[i+1] pairs with hist_cg_rnorm[i].
            hist_obj_pull, hist_con_pull
                objective-vs-constraint balance: norm of the linearized-objective RHS
                pull and of the total constraint pull ||Σ Aᵀy||.
            hist_con_pull_op
                list of per-operator constraint pulls ||Aᵀy||.
        """

    def get_warmstart(self) -> dict:
        """
        Capture a warm-start snapshot from the last solve().

        Returns the state at the returned (best-feasible) iterate, or at the final
        iterate if none was feasible, as a plain, picklable dict. Pass it to
        set_warmstart() before the next solve().

        Returns
        -------
        dict with keys:
            active : bool   -- False if no snapshot is available (no solve has run).
            N, Naxis : int
            dt : float
            X, fixer : 1-D numpy arrays (length N*Naxis): primal and free/fixed mask.
            ops : list of dicts, one per operator (ordered like get_op_names()):
                key : str             -- operator key (see get_op_keys()), used for matching.
                y : 1-D numpy array   -- dual (Lagrange multiplier), in normalized Ax-space.
                weight, gamma : float -- ADMM weight and relaxation.
                spec_norm : float     -- operator normalization at capture, used to rescale y.
                blocks : list[int]    -- Ax-space partition, used to resize y.
        """

    def set_warmstart(self, warmstart: dict) -> None:
        """
        Load a warm-start snapshot (from get_warmstart()) for the next solve().

        Applies to the next solve() only; later solves start cold from X0.

        Seeds the primal and each operator's dual, weight, and gamma. Across a grid
        change, free runs of the waveform and each dual block are resampled, and fixed
        runs come from the new set_vals. Operators are matched by key (see
        get_op_keys()), so add them in the same order; unmatched operators (e.g. a newly
        added constraint) start cold. If Naxis or the number of free runs per axis
        differs, the snapshot is ignored with a warning.

        Parameters
        ----------
        warmstart : dict
            A snapshot as returned by get_warmstart(). Requires keys N, Naxis, dt, X,
            fixer, and ops (each with key, y, weight, gamma, spec_norm, blocks).
        """

    def set_general_params(self, min_iter: int = 1, max_iter: int = 2000, log_interval: int = 20, gamma_x: float = 1.6, max_feval: int = 12000, obj_patience: int = 20) -> None:
        """
        Set general solver parameters.

        Deprecated: prefer setting the properties directly, e.g. `solver.max_iter = 5000`.

        Parameters
        ----------
        min_iter : int, optional
            Minimum outer iterations before feasible iterates are considered.
        max_iter : int, optional
            Maximum outer iterations.
        log_interval : int, optional
            Outer iterations between debug log prints.
        gamma_x : float, optional
            Over-relaxation factor for the outer primal update.
        max_feval : int, optional
            Maximum total inner (CG) iterations over the whole solve.
        obj_patience : int, optional
            Objective problems: stop after this many feasible iterations with no
            objective improvement and return the best feasible iterate. Without an
            objective the first feasible iterate is returned.
        """

    def set_ils_params(self, ils_tol: float = 0.001, ils_max_iter: int = 20, ils_min_iter: int = 2, ils_sigma: float = 0.0001, ils_tik_lam: float = 0.0) -> None:
        """
        Set inner (indirect) linear solver parameters.

        Deprecated: prefer setting the properties directly, e.g. `solver.ils_tol = 0.1`.

        Parameters
        ----------
        ils_tol : float, optional
            Inner-solver tolerance, relative to the warm-start residual.
        ils_max_iter : int, optional
            Maximum inner-solver iterations per outer iteration.
        ils_min_iter : int, optional
            Minimum inner-solver iterations per outer iteration.
        ils_sigma : float, optional
            Proximal weight: adds sigma*I to the inner system.
        ils_tik_lam : float, optional
            Tikhonov weight: adds tik_lam*I to the inner-system matrix.
        """

class SolverGroptSDMM(Solver):
    """SDMM solver for GrOpt gradient optimization problems."""

    def __init__(self) -> None: ...

    @property
    def bb_enable(self) -> bool:
        """
        Per-operator BB (spectral) reweighting every rw_interval iterations (default True).
        """

    @bb_enable.setter
    def bb_enable(self, arg: bool, /) -> None: ...

    @property
    def grw_enable(self) -> bool:
        """
        Global reweighting of the worst persistently infeasible operator (default True).
        """

    @grw_enable.setter
    def grw_enable(self, arg: bool, /) -> None: ...

    @property
    def rw_interval(self) -> int:
        """BB reweighting interval (iterations)."""

    @rw_interval.setter
    def rw_interval(self, arg: int, /) -> None: ...

    @property
    def rw_e_corr(self) -> float:
        """BB reweighting correlation threshold."""

    @rw_e_corr.setter
    def rw_e_corr(self, arg: float, /) -> None: ...

    @property
    def rw_eps(self) -> float:
        """BB reweighting numerical-stability epsilon."""

    @rw_eps.setter
    def rw_eps(self, arg: float, /) -> None: ...

    @property
    def rw_scalelim(self) -> float:
        """Max factor by which one BB update can raise or lower a weight."""

    @rw_scalelim.setter
    def rw_scalelim(self, arg: float, /) -> None: ...

    @property
    def grw_min_infeasible(self) -> int:
        """
        grw: consecutive infeasible iterations before an operator can be bumped.
        """

    @grw_min_infeasible.setter
    def grw_min_infeasible(self, arg: int, /) -> None: ...

    @property
    def grw_interval(self) -> int:
        """grw: reweighting interval (iterations)."""

    @grw_interval.setter
    def grw_interval(self, arg: int, /) -> None: ...

    @property
    def grw_mod(self) -> float:
        """grw: multiplicative weight-bump factor."""

    @grw_mod.setter
    def grw_mod(self, arg: float, /) -> None: ...

    @property
    def grw_balanced(self) -> bool:
        """
        grw: after a bump, divide all ADMM constraint weights by grw_mod^(1/K) (K = number of ADMM constraints) to keep their geometric mean fixed. Default False.
        """

    @grw_balanced.setter
    def grw_balanced(self, arg: bool, /) -> None: ...

    @property
    def reproject_iterate(self) -> bool:
        """
        Re-project the over-relaxed iterate onto the equality constraints each outer iteration (default True). NLCG/BiCGstabl require it when any constraint uses project=True.
        """

    @reproject_iterate.setter
    def reproject_iterate(self, arg: bool, /) -> None: ...

    @property
    def cutoff_freq(self) -> float:
        """
        Low-pass cutoff [Hz] applied to the iterate each outer iteration (per-free-run DST-I, see low_freq_project) to suppress high-frequency oscillation; <= 0 disables (default).
        """

    @cutoff_freq.setter
    def cutoff_freq(self, arg: float, /) -> None: ...

    @property
    def cutoff_iter(self) -> int:
        """
        Outer iteration at which the cutoff_freq projection stops; < 0 = project on every iteration.
        """

    @cutoff_iter.setter
    def cutoff_iter(self, arg: int, /) -> None: ...

    @property
    def cutoff_trans(self) -> float:
        """
        Raised-cosine roll-off width as a fraction of the cutoff bin; 0 = brick wall (default). See low_freq_project.
        """

    @cutoff_trans.setter
    def cutoff_trans(self, arg: float, /) -> None: ...

    @property
    def tr_enable(self) -> bool:
        """
        Trust-region step control (default False): re-solve a step with a larger proximal sigma when tr_monitor rejects it.
        """

    @tr_enable.setter
    def tr_enable(self, arg: bool, /) -> None: ...

    @property
    def tr_tol(self) -> float:
        """
        Trust-region reject threshold; <= 0 uses the monitor default (0.2 for 'linearization_error', 0.02 for 'feasibility', 0.5 for 'rel_step'). For 'linearization_error' it is the max relative SAFE model error ||true - linear|| / ||true||.
        """

    @tr_tol.setter
    def tr_tol(self, arg: float, /) -> None: ...

    @property
    def tr_bump(self) -> float:
        """Proximal-sigma multiplier applied on each rejected step (default 4)."""

    @tr_bump.setter
    def tr_bump(self, arg: float, /) -> None: ...

    @property
    def tr_max_reject(self) -> int:
        """
        Max re-solves per outer iteration before taking the most-damped step anyway (default 5).
        """

    @tr_max_reject.setter
    def tr_max_reject(self, arg: int, /) -> None: ...

    @property
    def tr_decay(self) -> float:
        """
        Sigma relaxation factor toward ils_sigma on an accepted step (default 0.5).
        """

    @tr_decay.setter
    def tr_decay(self, arg: float, /) -> None: ...

    @property
    def tr_monitor(self) -> str:
        """
        Signal that drives the trust region: 'linearization_error' (default; relative error of the frozen SAFE linearization), 'feasibility' (reject if the SAFE violation grows past max(previous, tr_tol)), or 'rel_step' (||dx||/||x||). Only SAFE operators report the first two signals. 'none' disables the trust region; any other value also disables it, with a warning.
        """

    @tr_monitor.setter
    def tr_monitor(self, arg: str, /) -> None: ...

    @property
    def obj_gate_enable(self) -> bool:
        """
        Feasibility-gated objective (default False): scale the objective pull by exp(-violation/obj_gate_scale) so the objective acts only near feasibility. Only SAFE (PNS/CNS) constraints report a violation.
        """

    @obj_gate_enable.setter
    def obj_gate_enable(self, arg: bool, /) -> None: ...

    @property
    def obj_gate_scale(self) -> float:
        """
        Violation scale of the objective gate, in SAFE units (fraction of the stimulation limit; default 0.05). Smaller keeps the objective gated until closer to feasibility.
        """

    @obj_gate_scale.setter
    def obj_gate_scale(self, arg: float, /) -> None: ...

    def set_sdmm_params(self, rw_interval: int = 8, rw_e_corr: float = 0.4, rw_eps: float = 1e-36, rw_scalelim: float = 1.5, grw_min_infeasible: int = 20, grw_interval: int = 20, grw_mod: float = 2.0) -> None:
        """
        Set SDMM reweighting parameters.

        Deprecated: prefer setting the properties directly, e.g. `solver.rw_interval = 16`.

        Parameters
        ----------
        rw_interval : int, optional
            BB reweighting interval (iterations).
        rw_e_corr : float, optional
            BB reweighting correlation threshold.
        rw_eps : float, optional
            BB reweighting numerical-stability epsilon.
        rw_scalelim : float, optional
            Max factor by which one BB update can raise or lower a weight.
        grw_min_infeasible : int, optional
            grw: consecutive infeasible iterations before an operator can be bumped.
        grw_interval : int, optional
            grw: reweighting interval (iterations).
        grw_mod : float, optional
            grw: multiplicative weight-bump factor.
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

def solve(params: GroptParams, min_iter: int = 1, max_iter: int = 2000, log_interval: int = 20, gamma_x: float = 1.6, max_feval: int = 12000, obj_patience: int = 20, ils_tol: float = 0.001, ils_max_iter: int = 20, ils_min_iter: int = 2, ils_sigma: float = 0.0001, ils_tik_lam: float = 0.0) -> SolveResult:
    """
    Solve a GrOpt problem with a new SolverGroptSDMM.

    For other options (reweighting, trust region, warm starts, debug histories),
    create a SolverGroptSDMM and set its properties directly.

    Parameters
    ----------
    params : GroptParams
        The problem definition.
    min_iter : int, optional
        Minimum outer iterations before feasible iterates are considered.
    max_iter : int, optional
        Maximum outer iterations.
    log_interval : int, optional
        Outer iterations between debug log prints.
    gamma_x : float, optional
        Over-relaxation factor for the outer primal update.
    max_feval : int, optional
        Maximum total inner (CG) iterations over the whole solve.
    obj_patience : int, optional
        Objective problems: stop after this many feasible iterations with no
        objective improvement.
    ils_tol : float, optional
        Inner-solver tolerance, relative to the warm-start residual.
    ils_max_iter : int, optional
        Maximum inner-solver iterations per outer iteration.
    ils_min_iter : int, optional
        Minimum inner-solver iterations per outer iteration.
    ils_sigma : float, optional
        Proximal weight: adds sigma*I to the inner system.
    ils_tik_lam : float, optional
        Tikhonov weight: adds tik_lam*I to the inner-system matrix.

    Returns
    -------
    SolveResult
        The optimization result.
    """

class NormType(enum.Enum):
    L2 = 0

    Inf = 1

class EqProjSolver(enum.Enum):
    LDLT = 0

    COD = 1

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
        Index of the operator in all_op (see GroptParams.get_op_names()).

    Returns
    -------
    float
        Estimated spectral norm.
    """

def get_SAFE(G: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')], dt: float, true_safe: bool = True, new_first_axis: int = 0, demo_params: bool = True, safe_params: object | None = None) -> Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')]:
    """
    Compute the SAFE (PNS) response for a single-axis gradient waveform.

    Parameters
    ----------
    G : np.ndarray
        Gradient waveform [T/m], single axis.
    dt : float
        Raster time in seconds.
    true_safe : bool, optional
        Currently ignored.
    new_first_axis : int, optional
        Use the SAFE parameters of this axis (0, 1, 2) for the waveform.
    demo_params : bool, optional
        Use the built-in demo SAFE parameters. Must be True if safe_params is None.
    safe_params : dict, optional
        Dictionary of SAFE parameters (see gropt.readasc).

    Returns
    -------
    np.ndarray
        SAFE response at each time point, as a fraction of the stimulation limit
        (1.0 = at the limit).
    """

def low_freq_project(x: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')], dt: float, cutoff_hz: float, fixer: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')] = ..., Naxis: int = 1, trans_frac: float = 0.0) -> Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')]:
    """
    Low-pass a waveform with a per-free-run DST-I projection.

    This is the filter SolverGroptSDMM applies when cutoff_freq > 0. Each maximal run
    of free samples is band-limited independently with a DST-I, which assumes the
    run is bounded by zeros; fixed samples are unchanged.

    Parameters
    ----------
    x : np.ndarray
        Waveform, length Naxis*N, axis-major.
    dt : float
        Raster time [s].
    cutoff_hz : float
        Cutoff frequency [Hz]; <= 0 returns x unchanged.
    fixer : np.ndarray, optional
        Free mask (1 = free, 0 = fixed), length Naxis*N. Empty (default) or any
        other length treats every sample as free.
    Naxis : int, optional
        Number of axes.
    trans_frac : float, optional
        Raised-cosine roll-off width as a fraction of the cutoff bin; 0 (default) =
        brick wall. A wider roll-off suppresses more near-cutoff oscillation but
        also strips harmonics that keep plateaus flat.

    Returns
    -------
    np.ndarray
        The projected copy of x.
    """

def test_eigen_assertions(test_type: int) -> None:
    """
    Trigger a deliberate Eigen error, to check whether Eigen assertions are compiled in.

    Assertions are enabled by building with GROPT_EIGEN_ASSERTIONS=ON; tests 1 and 3
    then trip an Eigen assertion, which aborts the process.

    Parameters
    ----------
    test_type : int
        1 = out-of-bounds read, 2 = raise RuntimeError unless new matrices are
        NaN-initialized (also set by GROPT_EIGEN_ASSERTIONS), 3 = size mismatch.
    """
