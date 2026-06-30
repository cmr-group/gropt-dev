from collections.abc import Callable
import enum
from typing import Annotated, overload

import numpy
from numpy.typing import NDArray


__build_date__: str = 'Jun 14 2026 15:21:56'

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

    @property
    def normalize_obj(self) -> bool:
        """
        If True, self-normalize the linearized objective direction so weight_mod sets a constant step magnitude (the pull no longer grows with ||AᵀA x||) — a forgiving rate knob for b-value maximization (endpoint is the constrained max regardless of value).
        """

    @normalize_obj.setter
    def normalize_obj(self, arg: bool, /) -> None: ...

    @property
    def eq_proj_solver(self) -> EqProjSolver:
        """
        Linear solver for the equality projection (project=True constraints). EqProjSolver.LDLT (default) is the original exact Gram solve -- fast/precise for well-conditioned sets (moments, concomitant) but blows up on singular/collinear rows. EqProjSolver.COD is a rank-revealing decomposition of Mhat that handles collinear rows (e.g. multiple eddy time-constants), keeps full precision (no Gram condition-squaring), and stays a true projection inside the CG. Switch to COD when LDLT can't handle your constraint set.
        """

    @eq_proj_solver.setter
    def eq_proj_solver(self, arg: EqProjSolver, /) -> None: ...

    @property
    def eq_proj_rcond(self) -> float:
        """
        Rank tolerance for eq_proj_solver = COD only (singular values below eq_proj_rcond * largest are dropped); ignored by LDLT. Larger drops more near-dependent rows; <= 0 uses Eigen's default threshold. Default 1e-10.
        """

    @eq_proj_rcond.setter
    def eq_proj_rcond(self, arg: float, /) -> None: ...

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

    def add_gmax_vec(self, gmax_vec: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')], rot_variant: bool = True, weight_mod: float = 1.0) -> None:
        """
        Add a per-location gradient amplitude constraint.

        Parameters
        ----------
        gmax_vec : np.ndarray
            Gradient amplitude limit at each location [T/m]. Length N is broadcast
            across all axes; length Naxis*N sets a per-axis limit. Only supported
            with rot_variant=True.
        rot_variant : bool, optional
            Must be True. A vector limit on the gradient magnitude across axes is
            not supported.
        weight_mod : float, optional
            Weighting factor for this constraint.
        """

    def add_concomitant(self, start_idx: int = 0, rot_variant: bool = True, weight_mod: float = 1.0, tol0: float = 0.1, target: float = 1.0, fix_gamma: bool = False, gamma_fix: float = 1.0, project: bool = False, as_objective: bool = False) -> None:
        """
        Add a concomitant constraint.

        Balances the gradient energy integral (sum of g^2 * dt) before and after the
        180: constrains the ratio pos/neg to `target` (default 1 = balanced) within tol0.
        This is a nonconvex quadratic constraint, so it is handled as a soft ADMM
        constraint (it cannot use the linear moment null-space projection).

        The prox is the true Euclidean projection onto the feasible band: if the ratio
        is inside the band it is left alone; otherwise the pre/post blocks are scaled by
        minimal displacement onto the nearest boundary cone (reduces to the arithmetic
        mean of the two block norms for exact equality).

        Parameters
        ----------
        start_idx : int, optional
            Starting index for the constraint (-1 = beginning).
        rot_variant : bool, optional
            If True, use rotationally invariant formulation.
        weight_mod : float, optional
            Weighting factor for this constraint.
        tol0 : float, optional
            Fractional tolerance on the pre/post energy ratio (feasible when
            |pos/neg - target| <= tol0). The prox aims for a cushioned (tighter) band so
            the solver reaches tol0 with margin. Note: the dt in the integral cancels
            in this ratio for uniform dt, so it does not change the numbers today --
            it keeps the energies physical and correct for non-uniform dt.
        target : float, optional
            Target pre/post energy ratio pos/neg (default 1.0 = balanced). Applies in all
            modes: the box/prox aims for pos/neg = target +/- tol, and the projection /
            objective forms enforce the linearized pos - target*neg = 0.
        fix_gamma : bool, optional
            If True, hold this operator's ADMM relaxation gamma fixed (the reweighter
            will not adapt it). The concomitant prox is nonconvex, so the reweighter's
            over-relaxation (gamma -> ~1.9) has no stability guarantee here; pin it to
            damp the resulting jitter.
        gamma_fix : float, optional
            The fixed gamma used when fix_gamma is True. 1.0 = vanilla ADMM (no
            over-relaxation); < 1.0 = under-relaxed / extra damping.
            Ignored when project=True (no ADMM path).
        project : bool, optional
            If True, enforce the constraint EXACTLY via the equality null-space
            projection instead of as a soft ADMM constraint. The nonconvex energy
            balance is linearized at the current iterate each outer iteration (SQP)
            and projected jointly with any projected moments, so it no longer uses
            the prox/consensus and the jitter from fighting the objective is removed.
            weight_mod/tol0/fix_gamma are then irrelevant (tol0 still sets the
            feasibility report band).
        as_objective : bool, optional
            If True, drive the balance from the objective (all_obj) path via an
            augmented Lagrangian: c(x)=pos-neg is linearized each outer iteration and
            enforced with a scalar dual (method of multipliers) that injects a rank-1
            PSD curvature into the CG LHS and a pull into the RHS (weight_mod = the
            penalty rho). It contributes to the solve but is NOT the optimization
            target (does not affect best-feasible scoring) and is NOT feasibility-
            gated -- inspect pos/neg on the result. Mutually exclusive with project
            (project is ignored when as_objective=True). The three modes mirror
            b-value: default = soft ADMM constraint, project = exact projector,
            as_objective = augmented-Lagrangian objective.
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

    def add_smax_vec(self, smax_vec: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')], rot_variant: bool = True, weight_mod: float = 1.0) -> None:
        """
        Add a per-location slew rate constraint.

        Parameters
        ----------
        smax_vec : np.ndarray
            Slew-rate limit at each location [T/m/s]. Length (N-1) is broadcast
            across all axes; length Naxis*(N-1) sets a per-axis limit. Only
            supported with rot_variant=True.
        rot_variant : bool, optional
            Must be True. A vector limit on the slew magnitude across axes is
            not supported.
        weight_mod : float, optional
            Weighting factor for this constraint.
        """

    def add_moment(self, order: float = 0, target: float = 0.0, tol: float = 1e-06, units: str = 'mT*ms/m', axis: int = 0, start_idx: int = -1, stop_idx: int = -1, ref_idx: int = 0, weight_mod: float = 1.0, project: bool = False) -> None:
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

    def add_eddy(self, lam: object, tol: float = 0.0001, weight_mod: float = 1.0, project: bool = False) -> None:
        """
        Add an eddy current constraint (residual eddy current at the end of the waveform).

        The target is always 0. With project=False (default) it is a soft box constraint
        |eddy| <= tol per time constant; with project=True the eddy current is driven to
        EXACTLY 0 via the null-space equality projection (like add_moment), enforced
        inside the CG step rather than as an ADMM penalty.

        Parameters
        ----------
        lam : float or np.ndarray
            Time constant(s) for eddy currents [seconds]. A single float is treated
            as a one-element array. One equality row is added per (axis, time constant).
        tol : float, optional
            Box half-width for the constraint form (ignored when project=True).
        weight_mod : float, optional
            Weighting factor for this constraint (constraint form only).
        project : bool, optional
            If True, enforce eddy == 0 exactly by null-space projection instead of the
            soft box constraint. Use a small number of time constants -- many exact
            eddy equalities (plus moments/concomitant) can over-constrain the free DOFs.
        """

    def add_bvalue(self, target: float = 100.0, tol: float = 1.0, start_idx0: int = -1, stop_idx0: int = -1, weight_mod: float = 1.0, mode: object = 2, max_scale: float = 1.01, as_objective: bool = False, linearize: bool = True) -> None:
        """
        Add a b-value term (constraint by default, or a maximization objective).

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
        as_objective : bool, optional
            If True, add b-value as a maximization objective (all_obj,
            obj_weight=-1) instead of a constraint.
        linearize : bool, optional
            Objective only. If True (default), the objective enters as a
            linearized gradient in the RHS (convex-concave/DCA), keeping the CG
            matrix positive-definite (recommended). If False, it enters as
            curvature in the LHS, which can make the system indefinite — for
            comparison only.
        """

    def add_diff_basin(self, window_time: float, eps_factor: float, gmax: float, weight_mod: float = 1.0, same_sign: bool = False) -> None:
        """
        Add a diffusion basin-orientation constraint (single 180 in the middle).

        Forces the gradient to FLIP sign across the 180 -- the structure of the global-max
        diffusion waveform -- and breaks the global +/- symmetry, so b-value maximization
        commits to the correct basin in a SINGLE solve. Two one-sided constraints on the
        mean gradient in a short window on each side of the 180 (located internally from
        the inv_vec sign flip; the fixed 180 RF block is skipped):

            mean(g, pre window)  >= +eps
            mean(g, post window) <= -eps,    eps = eps_factor * gmax

        Keep eps_factor small so the optimum (which has a strong flip) satisfies it with
        slack -- the constraint is then inactive at the optimum and does not distort the
        result; it only fences off the wrong basin during the descent. Soft ADMM
        constraint, not an equality projection.

        Parameters
        ----------
        window_time : float
            Window width on each side of the 180 [s] (e.g. 1e-3).
        eps_factor : float
            Min |mean window gradient| as a fraction of gmax (e.g. 0.07). Big enough to
            clear the ambiguous flat zone, small enough not to bias the lobe magnitude.
        gmax : float
            Gradient limit [T/m], used to scale eps = eps_factor * gmax.
        weight_mod : float, optional
            Weighting factor for this constraint.
        same_sign : bool, optional
            False (default): force OPPOSITE signs across the 180 (pre >= +eps, post <= -eps)
            -- the flip / global-max basin. True: force the SAME sign (pre >= +eps,
            post >= +eps) -- the no-flip / local basin. (Negative eps does NOT do this; it
            only loosens the bound, so this is a flag.)

        Notes
        -----
        Assumes a single axis / single 180 (the diffusion case). Default sign convention is
        pre >= +eps, post <= -eps; keep any X0 seed orientation consistent with that.
        """

    def add_TV(self, tv_lam: float = 0.0, weight_mod: float = 1.0, order: int = 1) -> None:
        """
        Add total-variation (L1) regularization on a finite difference of the gradient.

        This is a PENALTY (prox of the L1 norm), not a budget constraint: it minimizes
        tv_lam * ||difference||_1 rather than enforcing TV <= tol, and it does not gate
        feasibility. tv_lam is a pure regularization WEIGHT on the physical slew
        (order 1) / jerk (order 2): the prox divides by the operator's ADMM weight, so
        the effective penalty is exactly tv_lam regardless of rho (the reweighter can
        adapt rho for convergence without changing the result), and you do not have to
        re-tune it when weight_mod changes. It is a weight, not a hard threshold -- use a
        separate operator if you need a hard slew/jerk limit.

        Parameters
        ----------
        tv_lam : float, optional
            L1 regularization WEIGHT on the physical difference (must be > 0 to have
            effect). Larger = smoother.
        weight_mod : float, optional
            Seeds this operator's initial ADMM weight (rho). Convergence/balance knob
            only -- it no longer changes the penalty strength (use tv_lam for that).
        order : int, optional
            1 = TV of the gradient (||slew||_1): penalizes slew magnitude -> sparse
                slew / blocky gradients. Note this fights bang-bang ramps.
            2 = TV of the slew (||jerk||_1): penalizes CHANGES in slew -> piecewise-
                constant slew (clean bang-bang). This is the one that removes slew
                jitter without fighting the ramps; use a small tv_lam as a tie-breaker.
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

    def get_op_names(self) -> list[str]:
        """
        Return the constraint-operator names, in all_op order.

        The order matches the per-operator index in the solver debug
        histories (hist_weight, hist_gamma, hist_r_prim, hist_r_dual),
        so the returned list can be used directly as plot labels.
        """

class Solver:
    @property
    def extra_debug(self) -> bool:
        """
        If True, populate debug_solver with per-iteration history after solve().
        """

    @extra_debug.setter
    def extra_debug(self, arg: bool, /) -> None: ...

    @property
    def min_iter(self) -> int:
        """Minimum outer iterations before stopping."""

    @min_iter.setter
    def min_iter(self, arg: int, /) -> None: ...

    @property
    def max_iter(self) -> int:
        """Maximum outer iterations."""

    @max_iter.setter
    def max_iter(self, arg: int, /) -> None: ...

    @property
    def log_interval(self) -> int:
        """Iterations between debug log prints."""

    @log_interval.setter
    def log_interval(self, arg: int, /) -> None: ...

    @property
    def gamma_x(self) -> float:
        """Outer relaxation / over-relaxation factor."""

    @gamma_x.setter
    def gamma_x(self, arg: float, /) -> None: ...

    @property
    def max_feval(self) -> int:
        """Maximum total inner (CG) iterations."""

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
        Inner-solver relative tolerance (warm-start-relative; ~0.1 for inexact, tight for near-exact).
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
        """Proximal term added to the inner system (sigma I)."""

    @ils_sigma.setter
    def ils_sigma(self, arg: float, /) -> None: ...

    @property
    def ils_tik_lam(self) -> float:
        """Extra Tikhonov regularization for the inner system."""

    @ils_tik_lam.setter
    def ils_tik_lam(self, arg: float, /) -> None: ...

    def get_debug(self) -> dict:
        """
        Return debug history as a dict of lists.

        Only populated when extra_debug=True before solve(). Per-operator
        lists are ordered to match GroptParams.get_op_names().

        Returns
        -------
        dict with keys:
            hist_X, hist_Ax, hist_z, hist_y, hist_Aty
                list of 1-D numpy arrays, one per logged iteration.
            hist_weight, hist_gamma, hist_r_prim, hist_r_dual, hist_r_feas, hist_feas
                list (per iteration) of lists (per operator).
                hist_r_prim/hist_r_dual : Boyd ADMM primal/dual residuals.
                hist_r_feas/hist_feas   : relative infeasibility and binary feasible flag.
            hist_all_feas
                list of ints (per iteration): 1 if ALL operators were feasible that
                iteration, else 0. Plot against hist_bvalue to tell whether the
                objective is still climbing while feasibility flickers off (the
                returned best-feasible iterate is earlier) or has plateaued while
                feasible.
            hist_gamma_x
                list of floats, one per logged iteration.
            hist_bvalue
                achieved b-value per iteration (empty list if no b-value operator).
            best_feasible_iter
                int: the outer iteration whose iterate was actually returned (the
                best feasible one), or -1 if no feasible iterate was found. Use to
                mark the returned solution, e.g. hist_bvalue[best_feasible_iter].
            hist_cg_iter, hist_cg_rnorm0, hist_cg_rnorm, hist_cg_bnorm0
                inner linear-solver diagnostics, one entry per outer iteration:
                iterations taken, initial/final residual norm, and ||b||.
                (Generic to the ILS base class, so valid for CG/BiCGstabl/NLCG.)
            hist_cg_neg_curv
                CG only: 1 if a non-positive curvature direction (pAp<=0) was hit
                that outer iteration (indefinite system — watch this when running
                b-value as an objective with an unbalanced obj_weight).
            hist_cg_min_curv, hist_cg_max_curv
                CG only: min/max Rayleigh quotient (curvature) over the solve.
                hist_cg_min_curv -> 0 warns of near-singularity / blowup before
                hist_cg_neg_curv ever trips; max/min is a rough condition estimate.
                Only populated when extra_debug=True (the curvature dot is gated).
            hist_obj_pull, hist_con_pull
                objective-vs-constraint balance per iteration: ||g_obj|| (DCA
                objective RHS pull) and ||Σ Aᵀy|| (total constraint pull). The
                balance ratio is hist_obj_pull / hist_con_pull.
            hist_con_pull_op
                per-iteration list of per-operator ||Aᵀy|| (constraint pulls),
                ordered like get_op_names().
        """

    def get_warmstart(self) -> dict:
        """
        Capture a warm-start snapshot from the last solve.

        Returns the state at the returned (best-feasible) iterate as a plain, pickle-
        friendly dict (e.g. to ship to loky workers). Pass it to set_warmstart() before
        the next solve(). The consensus z is regenerated as z = A x on load, so it is
        not stored.

        Returns
        -------
        dict with keys:
            active : bool   -- False if no solve has run.
            N, Naxis : int
            dt : float
            X, fixer : 1-D numpy arrays (length N*Naxis): primal and source free/fixed mask.
            ops : list of dicts, one per operator (ordered like get_op_names()):
                key : str            -- operator unique_name, used for matching.
                y : 1-D numpy array  -- dual / Lagrange multiplier.
                weight, gamma : float
                blocks : list[int]   -- Ax-space partition (used to resize y).
        """

    def set_warmstart(self, warmstart: dict) -> None:
        """
        Load a warm-start snapshot (from get_warmstart()) for the next solve().

        One-shot: it applies to the next solve() only and is then cleared, so a later
        solve() on the same solver is cold (uses setvec_X0 / X0) unless you call this
        again. Chained sweeps still work since each step calls set_warmstart explicitly.

        The next solve() seeds the primal, per-operator duals, and weights from it,
        resizing across a grid change (free segments interpolate; each operator's Ax
        blocks interpolate). Operators are matched by unique_name; an operator with no
        match (e.g. a newly added constraint) starts cold. Build the new GroptParams in
        the same operator order so the auto-assigned names line up.

        Parameters
        ----------
        warmstart : dict
            A snapshot as returned by get_warmstart().
        """

    def set_general_params(self, min_iter: int = 1, max_iter: int = 2000, log_interval: int = 20, gamma_x: float = 1.6, max_feval: int = 12000, obj_patience: int = 20) -> None:
        """
        [DEPRECATED -- prefer setting the properties directly, e.g. solver.max_iter = 5000] Set general solver parameters.

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
        obj_patience : int, optional
            For objective problems, stop after this many feasible iterations with no
            objective improvement (returns the best feasible iterate). Ignored when
            there is no objective (returns on the first feasible iterate).
        """

    def set_ils_params(self, ils_tol: float = 0.001, ils_max_iter: int = 20, ils_min_iter: int = 2, ils_sigma: float = 0.0001, ils_tik_lam: float = 0.0) -> None:
        """
        [DEPRECATED -- prefer setting the properties directly, e.g. solver.ils_tol = 0.1] Set indirect linear solver parameters.

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
        """BB reweighting per-step scale limit."""

    @rw_scalelim.setter
    def rw_scalelim(self, arg: float, /) -> None: ...

    @property
    def grw_min_infeasible(self) -> int:
        """grw: minimum infeasible streak before adaptive reweighting kicks in."""

    @grw_min_infeasible.setter
    def grw_min_infeasible(self, arg: int, /) -> None: ...

    @property
    def grw_interval(self) -> int:
        """grw: adaptive reweighting interval."""

    @grw_interval.setter
    def grw_interval(self, arg: int, /) -> None: ...

    @property
    def grw_mod(self) -> float:
        """grw: multiplicative weight-bump factor."""

    @grw_mod.setter
    def grw_mod(self, arg: float, /) -> None: ...

    def set_sdmm_params(self, rw_interval: int = 8, rw_e_corr: float = 0.4, rw_eps: float = 1e-36, rw_scalelim: float = 1.5, grw_min_infeasible: int = 20, grw_interval: int = 20, grw_mod: float = 2.0) -> None:
        """
        [DEPRECATED -- prefer setting the properties directly, e.g. solver.rw_interval = 16] Set SDMM-specific parameters.

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

def solve(params: GroptParams, min_iter: int = 1, max_iter: int = 2000, log_interval: int = 20, gamma_x: float = 1.6, max_feval: int = 12000, obj_patience: int = 20, ils_tol: float = 0.001, ils_max_iter: int = 20, ils_min_iter: int = 2, ils_sigma: float = 0.0001, ils_tik_lam: float = 0.0) -> SolveResult:
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

def resize_warmstart(gparams: GroptParams, warmstart: dict) -> dict:
    """
    Resize a warm-start snapshot onto a target problem WITHOUT solving.

    Runs only the resizer used inside solve(): segment-aware interpolation of the
    primal X (free segments interpolate; fixed segments come from the target's
    set_vals) and block-wise interpolation of each operator's dual y. Use it to
    inspect the pre-optimization warm start; no equality projection or optimization
    is applied.

    Parameters
    ----------
    gparams : GroptParams
        The TARGET problem (built at the new size). Prepared internally if needed.
    warmstart : dict
        A snapshot from Solver.get_warmstart().

    Returns
    -------
    dict with keys:
        compatible : bool   -- False if Naxis or free-segment topology differ (X then returned as-is).
        X : 1-D numpy array -- the resized primal on the target grid.
        ops : list of dicts (one per target operator, ordered like get_op_names()):
            key : str, matched : bool, y : 1-D numpy array (resized dual; zeros if unmatched).
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
