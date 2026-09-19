"""GrOpt diffusion helpers.

Simplest call (gmax, smax, moment nulling and b-value only):

    from gropt.diffusion import DiffParams, solve
    res = solve(DiffParams(TE=60e-3, T_90=3e-3, T_180=5e-3, T_readout=16e-3))
    G = res["X"]
"""

from __future__ import annotations

import copy
import itertools
from dataclasses import dataclass, replace
from timeit import default_timer as timer

import numpy as np

import gropt
from gropt.readasc import asc_to_safe, get_random_safe_params


# ===========================================================================
# Problem definition
# ===========================================================================
@dataclass(frozen=True)
class DiffParams:
    """Diffusion problem definition."""

    # --- required timing [s] ---
    TE: float
    T_90: float
    T_180: float
    T_readout: float
    dt: float = 400e-6

    T_pre: float | None = None   # preencode only

    # --- gradient amplitude ---
    gmax: float = 0.08           # [T/m]
    w_gmax: float = 1.0

    # --- slew rate ---
    smax: float = 200.0          # [T/m/s]
    w_smax: float = 1.0

    # --- moments ---
    MMT: int = 0                 # null moments M0..M_MMT
    moment_project: bool = True  # exact null-space projection (recommended)
    moment_tol: float = 1e-5     # M0 tolerance; order k uses moment_tol * ||A_k|| / ||A_0||
    w_moment: float = 1.0

    # --- diffusion ---
    bvalue: float = 1000.0          # target b [s/mm^2] (te_search default)
    diff_mode: str = "gropt"        # "gropt" | "conventional" | "preencode"
    bval_mode: str = "obj"          # "obj" (maximize b) | constraint: "setval" | "minval" | "minval_max"
    bval_min: float = 100.0         # constraint modes: target b [s/mm^2]
    bval_obj_weight: float = 1.0    # obj mode: magnitude of the normalized b-value pull
    bval_max_scale: float = 1.02    # minval_max: per-iteration b scale factor
    w_bval: float = 1.0             # constraint modes only

    # --- SAFE, PNS, CNS ---
    pns_lim: float | None = None         # SAFE PNS limit (1.0 = model stim_limit); None = off
    cns_lim: float | None = None         # SAFE cardiac limit, same scale; None = off
    safe_params: SafeSource | None = None   # SAFE model source; None = SafeSource() (random, seed 42)
    safe_eps: float = 0.0                # softabs smoothing of SAFE |.| [T/m/s]; 0 = exact |.|
    w_pns: float = 1.0
    w_cns: float = 1.0

    # --- optional constraints (None / False => off) ---
    concomitant: bool = False
    w_concomitant: float = 1.0
    concomitant_project: bool = True

    eddy_lam: float | None = None   # eddy time constant [s]
    w_eddy: float = 1.0
    eddy_project: bool = True

    # --- basin control (None => off) ---
    basin_same_sign: bool | None = None  # False: force sign-flip (global) basin; True: no-flip (local)
    basin_window: float = 1e-3           # averaging window on each side of the 180 [s]
    basin_eps: float = 0.07              # min |mean g| per window, as a fraction of gmax

    # --- order-2 TV (jerk) regularizer: reduces slew jitter on the flat b-value optimum ---
    jerk_lam: float = 0.0           # add_TV(order=2) lambda, physical units (~1e-9..1e-6 useful); 0 = off
    w_jerk: float = 1.0

    # --- initial waveform (solve seed) ---
    x0_mode: str = "diff_init"      # "diff_init" (built-in 1e-2 seed) | "const" | "sine"
    x0_amp: float = 0.01            # seed amplitude on free samples (const, sine)
    x0_invert: bool = True          # flip the seed sign after the 180
    x0_periods: float = 1.0         # sine mode: periods per free run (zero at both ends)
    x0_project: bool = False        # pre-project the seed onto the moment null-space (M0..M_MMT = 0)


# ===========================================================================
# Solver settings
# ===========================================================================
@dataclass(frozen=True)
class SolverCfg:
    """Solver settings, tuned for diffusion b-value maximization.

    These intentionally differ from the general C++ / ``gropt.solve()`` defaults: inexact CG
    (``ils_tol=0.1``), a larger iteration budget, and more aggressive BB reweighting
    (``rw_e_corr=0.2``, ``rw_scalelim=2.0``).
    """

    # outer loop
    max_iter: int = 4000
    max_feval: int = 200000         # cap on total inner CG iterations
    min_iter: int = 1
    obj_patience: int = 20          # obj mode: stop after N feasible iters with no improvement
    obj_rtol: float = 1e-4          # relative improvement threshold for obj_patience
    gamma_x: float = 1.6            # over-relaxation of the X update

    # inner (CG) solver
    ils_tol: float = 0.1
    ils_max_iter: int = 20
    ils_min_iter: int = 2
    ils_sigma: float = 1e-4         # proximal weight
    ils_tik_lam: float = 0.0

    bb_reweight: bool = True        # BB per-operator adaptation
    rw_interval: int = 8            # BB cadence: reweight every rw_interval outer iters
    rw_e_corr: float = 0.2          # BB correlation threshold (Xu et al.); lower = more frequent updates
    rw_scalelim: float = 2.0        # max weight change per BB step (1.0 = frozen)
    rw_eps: float = 1e-36           # safeguard floor for the correlation/norm guards

    grw: bool = True                # global reweighting: bump the most persistently infeasible op
    grw_interval: int = 20          # grw cadence (only when grw)
    grw_mod: float = 2.0            # grw multiplier (only when grw)
    grw_balanced: bool = False      # rescale all weights to keep their geometric mean fixed

    reproject_iterate: bool = True  # re-project onto the equality constraints each outer iteration

    # low-pass the iterate each outer iteration (suppresses high-frequency oscillation)
    cutoff_freq: float = -1.0       # [Hz]; <= 0 off
    cutoff_iter: int = -1           # apply until this iteration; -1 = always
    cutoff_trans: float = 0.0       # raised-cosine roll-off, fraction of the cutoff bin; 0 = brick wall

    # trust region: re-solve with a larger proximal sigma when the step monitor rejects a step
    tr_enable: bool = False
    tr_tol: float = -1.0            # reject threshold; <= 0 uses the monitor default
    tr_bump: float = 4.0            # sigma multiplier per reject
    tr_max_reject: int = 5          # re-solves per outer iteration before taking the step anyway
    tr_decay: float = 0.5           # sigma decay toward ils_sigma on an accepted step
    tr_monitor: str = "linearization_error"  # | "feasibility" | "rel_step"

    # feasibility-gated objective: b-pull *= exp(-violation/obj_gate_scale); only SAFE reports violation
    obj_gate: bool = False
    obj_gate_scale: float = 0.05

    extra_debug: bool = False       # populate solver.get_debug()


# ===========================================================================
# SAFE (PNS / cardiac) parameter source
# ===========================================================================
@dataclass(frozen=True)
class SafeSource:
    """Hashable pointer to SAFE model parameters, resolved to ``(pns, cns)`` dicts by ``resolve``.

    * ``kind="random"`` -- synthetic params from ``seed`` (default).
    * ``kind="asc"``    -- params parsed from the scanner file ``asc_file`` (required).
    """

    kind: str = "random"
    seed: int = 42
    asc_file: str | None = None

    def resolve(self):
        """Resolve the pointer to concrete SAFE parameter dicts.

        Returns
        -------
        tuple of dict
            ``(params_pns, params_cns)`` -- the PNS and cardiac SAFE model parameters.

        Raises
        ------
        ValueError
            If ``kind='asc'`` without ``asc_file``, or ``kind`` is unrecognized.
        """
        if self.kind == "random":
            return get_random_safe_params(self.seed)
        if self.kind == "asc":
            if self.asc_file is None:
                msg = "SafeSource(kind='asc') requires asc_file"
                raise ValueError(msg)
            return asc_to_safe(self.asc_file)
        msg = f"Unknown SafeSource kind: {self.kind!r}"
        raise ValueError(msg)



# ===========================================================================
# TE search
# ===========================================================================
def te_search(base_cfg: DiffParams, scfg: SolverCfg = None, target_b: float | None = None,
              mode: str = "minval", te_lo: float | None = None, te_hi: float | None = None,
              te_tol: float | None = None, max_expand: int = 6, verbose: bool = False):
    """Smallest TE whose waveform reaches ``b >= target_b``, by monotonic bisection.

    Modes:

    * ``"minval"`` (default): b-value constraint ``bval_min=target_b``; feasible if the solve converges.
      Robust under tight constraints, but infeasible TEs run to ``max_iter`` and are slow.
    * ``"obj"``: maximize b; feasible if converged and ``bvalue >= target_b``. Fast everywhere, but
      b-maximization can fail under tight constraints (e.g. strict PNS), which over-estimates TE.

    Parameters
    ----------
    base_cfg : DiffParams
        Problem definition. Its ``TE`` is swept; ``bval_mode``/``bval_min`` are overridden per ``mode``.
    scfg : SolverCfg, optional
        Solver settings for each evaluation (defaults to ``SolverCfg()``).
    target_b : float, optional
        b-value to reach [s/mm^2]; defaults to ``base_cfg.bvalue``.
    mode : str, optional
        ``"minval"`` (default) or ``"obj"``, see above.
    te_lo, te_hi : float, optional
        Initial TE bracket [s]. Defaults: ``te_lo`` = the minimum physical TE (b~0 there, infeasible),
        ``te_hi`` = minimum physical TE + 80 ms.
    te_tol : float, optional
        TE resolution to stop at [s]. Defaults to ``dt`` and is clamped to at least ``dt``, since TE is
        quantized to the dt grid.
    max_expand : int, optional
        Maximum number of times the bracket span is doubled while ``te_hi`` is infeasible. Default 6.
    verbose : bool, optional
        Print each evaluation.

    Returns
    -------
    dict with keys:
        TE        : smallest feasible TE [s] (bin start; None if none found in the expanded bracket).
        N         : the sample count at that TE (``int((TE - T_readout)/dt) + 1``), or None.
        cfg       : ``base_cfg`` at the found TE (original bval_mode etc.), e.g. for
                    ``plot_diff(out["cfg"], out["result"])`` or a re-solve. None if not found.
        result    : the solve dict at that TE (the feasibility waveform), or None.
        n_solves  : number of solves (evaluations in the same dt bin are cached).
        bracket   : (te_lo, te_hi) final bracket.
        search_time : wall-clock seconds for the whole search.
    """
    t_start = timer()
    scfg = scfg or SolverCfg()
    target_b = base_cfg.bvalue if target_b is None else target_b
    if mode not in ("minval", "obj"):
        msg = f"te_search mode must be 'minval' or 'obj', got {mode!r}"
        raise ValueError(msg)
    dt = base_cfg.dt
    t_ro = base_cfg.T_readout

    # TE is quantized to the dt grid: N = int((TE - T_readout)/dt) + 1
    te_tol = dt if te_tol is None else max(te_tol, dt)

    def grid_index(te):
        return int((te - t_ro) / dt)          # matches diff_init's int() truncation

    def bin_start(m):
        return t_ro + m * dt                  # smallest TE giving grid index m (N = m + 1)

    def bin_solve_te(m):
        return t_ro + (m + 0.5) * dt          # canonical in-bin TE (robustly maps back to m)

    cache = {}  # grid index m -> (feasible?, result)
    n_solves = 0

    def feasible(te):
        nonlocal n_solves
        m = grid_index(te)
        if m not in cache:
            n_solves += 1
            if mode == "minval":
                # b >= target_b as a constraint; a converged solve is feasible
                cfg = replace(base_cfg, TE=bin_solve_te(m), bval_mode="minval", bval_min=target_b)
            else:
                # maximize b; a blow-up leaves converged=False and reads as infeasible
                cfg = replace(base_cfg, TE=bin_solve_te(m), bval_mode="obj")
            res = solve(cfg, scfg)
            ok = bool(res["converged"]) and res["bvalue"] >= 0.999 * target_b
            cache[m] = (ok, res)
            if verbose:
                print(f"  TE={bin_start(m) * 1e3:6.2f} ms (N={m + 1}) -> b={res['bvalue']:8.1f} "
                      f"conv={res['converged']!s:5} {'feasible' if ok else 'infeasible'}", flush=True)
        return cache[m]

    # physical minimum TE (b ~ 0 there -> infeasible lower end)
    min_te = max(base_cfg.T_90 + base_cfg.T_180 + base_cfg.T_readout,
                 2.0 * (base_cfg.T_180 + base_cfg.T_readout))
    if te_lo is None:
        te_lo = min_te
    if te_hi is None:
        te_hi = min_te + 80e-3

    def out(m, res):
        # cfg keeps base_cfg's bval_mode etc., at the bin-start TE
        return {"TE": bin_start(m), "N": m + 1, "cfg": replace(base_cfg, TE=bin_start(m)),
                "result": res, "n_solves": n_solves, "bracket": (te_lo, te_hi),
                "search_time": timer() - t_start}

    # ensure the high end is feasible, growing the span if not
    ok_hi, res_hi = feasible(te_hi)
    for _ in range(max_expand):
        if ok_hi:
            break
        te_hi = te_lo + 2.0 * (te_hi - te_lo)
        ok_hi, res_hi = feasible(te_hi)
    if not ok_hi:
        return {"TE": None, "N": None, "cfg": None, "result": None, "n_solves": n_solves,
                "bracket": (te_lo, te_hi), "search_time": timer() - t_start}

    # low end already feasible: the answer is at/below it
    ok_lo, res_lo = feasible(te_lo)
    if ok_lo:
        return out(grid_index(te_lo), res_lo)

    # bisection: invariant te_lo infeasible, te_hi feasible; stop at the dt-grid resolution
    best_m, best_res = grid_index(te_hi), res_hi
    while te_hi - te_lo > te_tol:
        te_mid = 0.5 * (te_lo + te_hi)
        ok, res = feasible(te_mid)
        if ok:
            te_hi, best_m, best_res = te_mid, grid_index(te_mid), res
        else:
            te_lo = te_mid

    return out(best_m, best_res)


def solve(cfg: DiffParams, scfg: SolverCfg = None, warmstart: dict = None, keep_weights: bool = True,
          return_solver: bool = False):
    """Build + solve a single ``DiffParams``. The main entry point.

    Parameters
    ----------
    cfg : DiffParams
        Problem definition.
    scfg : SolverCfg, optional
        Solver settings; defaults to ``SolverCfg()`` (BB + global reweighting on, inexact CG).
    warmstart : dict, optional
        A snapshot from a previous result's ``"warmstart"`` to seed this solve.
    keep_weights : bool, optional
        Only relevant with a ``warmstart``. True (default) inherits each operator's ADMM weight from
        the snapshot; False resets them to this ``cfg``'s ``weight_mod`` (keeping the primal and
        duals). Use False for weight sweeps and continuation.
    return_solver : bool, optional
        If True, also return ``(solver, gp)`` for interactive inspection. The default keeps the return
        value picklable for parallel workers.

    Returns
    -------
    dict  (or ``(dict, solver, gp)`` if ``return_solver``)
        Keys: TE, dt, start_idx, bvalue, converged, n_iter, n_feval, X, warmstart, solve_time
        [, debug, op_names].
    """
    scfg = scfg or SolverCfg()
    gp, start_idx, op_weights = build_gparams(cfg)
    solver = make_solver(scfg)
    if warmstart is not None:
        ws = warmstart if keep_weights else apply_fresh_weights(warmstart, op_weights)
        solver.set_warmstart(ws)
    start_t = timer()
    r = solver.solve(gp)
    solve_t = timer() - start_t
    out = _result_dict(cfg, r, solver, scfg, gp, start_idx)
    out["solve_time"] = solve_t
    if return_solver:
        return out, solver, gp
    return out


# ===========================================================================
# Builders
# ===========================================================================
def _build_x0(gp, cfg: DiffParams, start_idx: int):
    """Construct an initial waveform per ``cfg.x0_mode``.

    Fixed regions (90/180 gaps, endpoints) always take their set_vals (0).

    Parameters
    ----------
    gp : gropt.GroptParams
        Prepared params providing the inv/fixer/set-value vectors and the grid (``N``, ``dt``).
    cfg : DiffParams
        Problem definition; the ``x0_*`` fields drive the seed.
    start_idx : int
        First free index (moment/inv_vec offset; ``> 0`` only for preencode).

    Returns
    -------
    numpy.ndarray or None
        The flat initial waveform, or None to keep the built-in ``diff_init`` seed.
    """
    if cfg.x0_mode == "diff_init":
        return None
    N, dt = gp.N, gp.dt
    inv = np.asarray(gp.getvec_inv_vec())
    fixer = np.asarray(gp.getvec_fixer())     # 1 = free, 0 = fixed
    setv = np.asarray(gp.getvec_set_vals())   # NaN = free, value = fixed
    free = fixer > 0.5
    x = np.zeros(N)

    if cfg.x0_mode == "const":
        x[free] = cfg.x0_amp
        if cfg.x0_invert:
            x = x * inv                        # +amp pre-180, -amp post-180
    elif cfg.x0_mode == "sine":
        # each free run gets x0_periods sine periods, zero at both ends (continuous with the fixed zeros)
        i = 0
        while i < N:
            if not free[i]:
                i += 1
                continue
            j = i
            while j < N and free[j]:
                j += 1
            n = j - i
            if n >= 2:
                s = np.arange(n) / (n - 1)
                x[i:j] = cfg.x0_amp * np.sin(2.0 * np.pi * cfg.x0_periods * s)
            i = j
        if cfg.x0_invert:
            x[inv < 0] *= -1.0                 # flip the post-180 lobes (post_sign = -1)
    else:
        msg = f"Unknown x0_mode: {cfg.x0_mode!r}"
        raise ValueError(msg)

    x[~free] = np.nan_to_num(setv[~free])      # fixed regions -> their set_vals

    if cfg.x0_project and cfg.MMT >= 0:
        # project onto the moment null-space (M0..M_MMT = 0), moving only free DOFs (like the solver)
        t = np.arange(N) * dt
        M = np.array([(t ** k) * inv for k in range(cfg.MMT + 1)])
        M[:, :start_idx] = 0.0
        Mt = M * fixer[None, :]                # free-DOF masked
        x = x - Mt.T @ np.linalg.lstsq(Mt @ Mt.T, M @ x, rcond=None)[0]
    return x


def suggest_moment_tol(cfg: DiffParams, cushion: float = 10.0):
    """Estimate the moment-nulling numerical floor and suggest a ``moment_tol``.

    Parameters
    ----------
    cfg : DiffParams
        Problem definition (uses ``dt``, ``gmax``, and the timing that sets the free window).
    cushion : float, optional
        Safety factor over the floor. Default 10 (one order of magnitude).

    Returns
    -------
    dict
        Keys: ``floor_M0``, ``suggested_tol``, ``T_free``, ``N_free``, ``A0_norm``, ``G_norm_est``.
    """
    gp = gropt.GroptParams()
    if cfg.diff_mode == "conventional":
        gp.diff_init_deadtime(dt=cfg.dt, TE=cfg.TE, T_90=cfg.T_90, T_180=cfg.T_180, T_readout=cfg.T_readout)
    elif cfg.diff_mode == "preencode":
        gp.diff_init_preencode(dt=cfg.dt, TE=cfg.TE, T_90=cfg.T_90, T_180=cfg.T_180,
                               T_readout=cfg.T_readout, T_pre=cfg.T_pre or 0.0)
    else:
        gp.diff_init(dt=cfg.dt, TE=cfg.TE, T_90=cfg.T_90, T_180=cfg.T_180, T_readout=cfg.T_readout)

    dt = gp.dt
    n_free = int(np.count_nonzero(np.asarray(gp.getvec_inv_vec())))
    t_free = dt * n_free
    a0_norm = 1e6 * dt * np.sqrt(n_free)          # ||A_0|| over the free window
    g_norm_est = cfg.gmax * np.sqrt(n_free)       # ||G|| ~ bang-bang near gmax
    floor = np.finfo(float).eps * a0_norm * g_norm_est   # = eps * 1e6 * gmax * T_free
    return {"floor_M0": floor, "suggested_tol": cushion * floor, "T_free": t_free,
                "N_free": n_free, "A0_norm": a0_norm, "G_norm_est": g_norm_est}


def build_gparams(cfg: DiffParams):
    """Build a ``GroptParams`` from a ``DiffParams``.

    Parameters
    ----------
    cfg : DiffParams
        Problem definition.

    Returns
    -------
    gp : gropt.GroptParams
        The prepared params object.
    start_idx : int
        First free index (0 for gropt/conventional; ``> 0`` for preencode).
    op_weights : dict
        Maps each operator's warm-start key to the ``weight_mod`` this cfg built it with
        (see :func:`apply_fresh_weights`).
    """
    gp = gropt.GroptParams()

    start_idx = 0
    if cfg.diff_mode == "gropt":
        gp.diff_init(dt=cfg.dt, TE=cfg.TE, T_90=cfg.T_90, T_180=cfg.T_180, T_readout=cfg.T_readout)
    elif cfg.diff_mode == "conventional":
        gp.diff_init_deadtime(dt=cfg.dt, TE=cfg.TE, T_90=cfg.T_90, T_180=cfg.T_180, T_readout=cfg.T_readout)
    elif cfg.diff_mode == "preencode":
        if cfg.T_pre is None:
            msg = "diff_mode='preencode' requires T_pre"
            raise ValueError(msg)
        start_idx = gp.diff_init_preencode(dt=cfg.dt, TE=cfg.TE, T_90=cfg.T_90, T_180=cfg.T_180,
                                           T_readout=cfg.T_readout, T_pre=cfg.T_pre)
    else:
        msg = f"Unknown diff_mode: {cfg.diff_mode!r}"
        raise ValueError(msg)

    # initial weight_mod per constraint operator, in all_op order (keyed after prepare())
    op_weights = []

    # always-on hardware + moment nulling
    gp.add_gmax(cfg.gmax, weight_mod=cfg.w_gmax)
    op_weights.append(cfg.w_gmax)

    gp.add_smax(cfg.smax, weight_mod=cfg.w_smax)
    op_weights.append(cfg.w_smax)

    for m in range(cfg.MMT + 1):
        gp.add_moment(m, 0.0, start_idx=start_idx, tol=cfg.moment_tol,
                      weight_mod=cfg.w_moment, project=cfg.moment_project)
        op_weights.append(cfg.w_moment)

    # SAFE (PNS / cardiac), both from one source (default: random params)
    if cfg.pns_lim is not None or cfg.cns_lim is not None:
        pns_params, cns_params = (cfg.safe_params or SafeSource()).resolve()
        gp.safe_eps = cfg.safe_eps  # copied into each Op_SAFE by add_SAFE; must be set first
        if cfg.pns_lim is not None:
            gp.add_SAFE(cfg.pns_lim, safe_params=pns_params, weight_mod=cfg.w_pns)
            op_weights.append(cfg.w_pns)
        if cfg.cns_lim is not None:
            gp.add_SAFE(cfg.cns_lim, safe_params=cns_params, weight_mod=cfg.w_cns)
            op_weights.append(cfg.w_cns)

    if cfg.concomitant:
        gp.add_concomitant(start_idx=start_idx, project=cfg.concomitant_project, weight_mod=cfg.w_concomitant)
        op_weights.append(cfg.w_concomitant)
    if cfg.eddy_lam is not None:
        gp.add_eddy(cfg.eddy_lam, weight_mod=cfg.w_eddy, project=cfg.eddy_project)
        op_weights.append(cfg.w_eddy)

    # order-2 TV (jerk) stability regularizer
    if cfg.jerk_lam > 0.0:
        gp.add_TV(cfg.jerk_lam, weight_mod=cfg.w_jerk, order=2)
        op_weights.append(cfg.w_jerk)

    # basin orientation (opt-in): keep small, it is inactive at the optimum
    if cfg.basin_same_sign is not None:
        gp.add_diff_basin(cfg.basin_window, cfg.basin_eps, cfg.gmax, same_sign=cfg.basin_same_sign)
        op_weights.append(1.0)  # add_diff_basin weight_mod default

    # b-value: objective (maximize) or constraint
    if cfg.bval_mode == "obj":
        gp.add_bvalue(as_objective=True, start_idx0=start_idx, weight_mod=cfg.bval_obj_weight)
        gp.normalize_obj = True
    else:
        gp.add_bvalue(cfg.bval_min, mode=cfg.bval_mode, start_idx0=start_idx,
                      weight_mod=cfg.w_bval, max_scale=cfg.bval_max_scale)
        op_weights.append(cfg.w_bval)

    # optional custom seed (None keeps the diff_init seed)
    x0 = _build_x0(gp, cfg, start_idx)
    if x0 is not None:
        gp.setvec_X0(x0, set_others=False)

    gp.prepare()

    # key the weights by warm-start key; the count check catches an add_* without its op_weights entry
    names = gp.get_op_keys()
    if len(names) != len(op_weights):
        msg = (f"tracked op_weights ({len(op_weights)}) != built operators ({len(names)}: "
               f"{', '.join(names)}); an add_* call in build_gparams is missing its op_weights entry")
        raise RuntimeError(msg)
    return gp, start_idx, dict(zip(names, op_weights, strict=True))


def make_solver(scfg: SolverCfg | None = None):
    """Build and configure a ``SolverGroptSDMM`` from a ``SolverCfg``.

    Parameters
    ----------
    scfg : SolverCfg, optional
        Solver settings; defaults to ``SolverCfg()``.

    Returns
    -------
    gropt.SolverGroptSDMM
        The configured solver.
    """
    if scfg is None:
        scfg = SolverCfg()

    s = gropt.SolverGroptSDMM()
    s.max_iter = scfg.max_iter
    s.min_iter = scfg.min_iter
    s.obj_patience = scfg.obj_patience
    s.obj_rtol = scfg.obj_rtol
    s.gamma_x = scfg.gamma_x
    s.max_feval = scfg.max_feval
    s.extra_debug = scfg.extra_debug
    s.cutoff_freq = scfg.cutoff_freq
    s.cutoff_iter = scfg.cutoff_iter
    s.cutoff_trans = scfg.cutoff_trans
    s.tr_enable = scfg.tr_enable
    s.tr_tol = scfg.tr_tol
    s.tr_bump = scfg.tr_bump
    s.tr_max_reject = scfg.tr_max_reject
    s.tr_decay = scfg.tr_decay
    s.tr_monitor = scfg.tr_monitor
    s.obj_gate_enable = scfg.obj_gate
    s.obj_gate_scale = scfg.obj_gate_scale
    s.reproject_iterate = scfg.reproject_iterate

    s.ils_tol = scfg.ils_tol
    s.ils_max_iter = scfg.ils_max_iter
    s.ils_min_iter = scfg.ils_min_iter
    s.ils_sigma = scfg.ils_sigma
    s.ils_tik_lam = scfg.ils_tik_lam

    s.bb_enable = scfg.bb_reweight
    s.rw_interval = scfg.rw_interval
    s.rw_e_corr = scfg.rw_e_corr
    s.rw_eps = scfg.rw_eps
    s.rw_scalelim = scfg.rw_scalelim
    s.grw_enable = scfg.grw
    s.grw_interval = scfg.grw_interval
    s.grw_mod = scfg.grw_mod
    s.grw_balanced = scfg.grw_balanced
    return s


# ===========================================================================
# Warm-start prototyping helpers
# ===========================================================================
def warmstart_summary(warmstart: dict):
    """Build a human-readable view of what a snapshot carries.

    Parameters
    ----------
    warmstart : dict
        A warm-start snapshot (as returned in a result's ``"warmstart"`` key).

    Returns
    -------
    list of dict
        One row per operator (in all_op order): ``key``, ``y_norm`` (dual magnitude that gets
        carried), ``weight`` (the ADMM rho the next solve would inherit), ``gamma``, and ``nblocks``.
    """
    rows = []
    for op in warmstart.get("ops", []):
        y = np.asarray(op["y"])
        rows.append({
            "key": op["key"],
            "y_norm": float(np.linalg.norm(y)),
            "weight": float(op["weight"]),
            "gamma": float(op["gamma"]),
            "nblocks": len(op.get("blocks", [])),
        })
    return rows


def apply_fresh_weights(warmstart: dict, op_weights: dict) -> dict:
    """Return a copy of ``warmstart`` with the carried weights replaced by freshly built ones.

    Operators are matched by ``key``, as ``set_warmstart`` does; snapshot operators missing from
    ``op_weights`` are left unchanged.

    Parameters
    ----------
    warmstart : dict
        The snapshot to copy and rewrite.
    op_weights : dict
        Maps each operator's warm-start key to its freshly built ``weight_mod``.

    Returns
    -------
    dict
        A deep copy of ``warmstart`` with each matched operator's ``weight`` replaced.
    """
    ws = copy.deepcopy(warmstart)
    for op in ws.get("ops", []):
        w = op_weights.get(op["key"])
        if w is not None:
            op["weight"] = float(w)
    return ws


# ===========================================================================
# Run helpers
# ===========================================================================
def _result_dict(cfg: DiffParams, r, solver, scfg: SolverCfg, gp=None, start_idx=0):
    """Assemble a plain, picklable result dict, always carrying a warm-start snapshot for chaining.

    Parameters
    ----------
    cfg : DiffParams
        Problem definition the result came from.
    r : object
        The raw solver result (``bvalue``, ``converged``, ``n_iter``, ``n_feval``, ``X``).
    solver : gropt.SolverGroptSDMM
        The solver, queried for the warm-start snapshot (and debug info if enabled).
    scfg : SolverCfg
        Solver settings; ``extra_debug`` gates the debug payload.
    gp : gropt.GroptParams, optional
        Params object, used only to attach ``op_names`` when debugging.
    start_idx : int, optional
        Moment/inv_vec offset carried into the result (0 for gropt/conventional).

    Returns
    -------
    dict
        Keys: TE, dt, start_idx, bvalue, converged, n_iter, n_feval, X, warmstart
        [, debug, op_names].
    """
    out = {
        "TE": cfg.TE,
        "dt": cfg.dt,
        "start_idx": start_idx,
        "bvalue": float(r.bvalue),
        "converged": bool(r.converged),
        "n_iter": int(r.n_iter),
        "n_feval": int(r.n_feval),
        "X": np.asarray(r.X),
        "warmstart": solver.get_warmstart(),
    }
    if scfg.extra_debug:
        out["debug"] = solver.get_debug()
        if gp is not None:
            out["op_names"] = gp.get_op_names()
    return out

def continuation(steps, base_cfg: DiffParams, scfg: SolverCfg = None, keep_weights: bool = False):
    """Run a sequential, warm-started continuation over a list of ``DiffParams`` overrides.

    Each step is a dict of overrides applied to ``base_cfg`` (e.g. a coarse-to-fine dt ramp
    ``[{"dt": 8e-4}, {"dt": 4e-4}, {"dt": 2e-4}]`` or a soft-to-stiff weight ramp) and warm-starts
    from the previous step's snapshot, resized across grid changes. Steps may change the operator
    set; a newly added operator starts cold.

    Parameters
    ----------
    steps : list of dict
        Per-step field overrides applied to ``base_cfg`` via ``replace``.
    base_cfg : DiffParams
        The template each step overrides.
    scfg : SolverCfg, optional
        Solver settings for every step (defaults to ``SolverCfg()``).
    keep_weights : bool, optional
        False (default) resets carried weights to each step's ``weight_mod``; True carries the adapted
        weights (useful for a pure dt ramp).

    Returns
    -------
    list of dict
        The per-step result dicts (each as from :func:`solve`).
    """
    scfg = scfg or SolverCfg()
    results = []
    ws = None
    for step in steps:
        res = solve(replace(base_cfg, **step), scfg, warmstart=ws, keep_weights=keep_weights)
        results.append(res)
        ws = res["warmstart"]
    return results


def dt_schedule(dt_target: float, n_levels: int = 3, factor: float = 2.0):
    """Build coarse->fine dt continuation steps: ``dt_target * factor**(n-1) ... dt_target``.

    Feed straight to :func:`continuation`, e.g.
    ``continuation(dt_schedule(2e-4, n_levels=3), cfg)``.

    Parameters
    ----------
    dt_target : float
        The finest (final) dt [s].
    n_levels : int, optional
        Number of levels, coarsest to finest. Default 3.
    factor : float, optional
        Per-level coarsening ratio. Default 2.0.

    Returns
    -------
    list of dict
        One ``{"dt": ...}`` override per level, coarsest first.
    """
    return [{"dt": dt_target * factor ** k} for k in range(n_levels - 1, -1, -1)]


# ===========================================================================
# Parallel sweeps
# ===========================================================================
def _grid_points(base_cfg, base_scfg, cfg_grid, scfg_grid):
    """Build the cross product of the cfg/scfg field-override grids.

    Parameters
    ----------
    base_cfg : DiffParams
        Template for the cfg axis.
    base_scfg : SolverCfg
        Template for the scfg axis.
    cfg_grid, scfg_grid : dict or None
        Field-override grids ``{field: [values]}``; ``None`` => that object is not varied.

    Returns
    -------
    list of tuple
        One ``(cfg, scfg)`` point per grid combination.
    """

    def combos(grid):
        if not grid:
            return [{}]
        keys = list(grid.keys())
        return [dict(zip(keys, vals, strict=True)) for vals in itertools.product(*(grid[k] for k in keys))]

    return [(replace(base_cfg, **c), replace(base_scfg, **s))
            for c in combos(cfg_grid) for s in combos(scfg_grid)]


def _reusable_executor(max_workers, pool_timeout):
    """Get loky's reusable pool with a long idle timeout (see ``pool_timeout`` in :func:`sweep`).

    Pass the same arguments each call so loky reuses the pool.

    Parameters
    ----------
    max_workers : int or None
        Worker count; None lets loky choose (~cpu count).
    pool_timeout : float
        Idle seconds before the pool is reaped.

    Returns
    -------
    loky reusable executor
        The (possibly pre-existing) reusable process pool.
    """
    from loky import get_reusable_executor

    return get_reusable_executor(max_workers=max_workers, timeout=pool_timeout)


def warm_pool(max_workers=None, pool_timeout=300.0):
    """Pre-spawn the loky worker pool so the one-time start-up cost is paid up front.

    No-op if loky is not installed. Use the same ``pool_timeout`` as :func:`sweep` so the pool is
    reused.

    Parameters
    ----------
    max_workers : int, optional
        Worker count; None lets loky choose (~cpu count).
    pool_timeout : float, optional
        Idle seconds before the pool is reaped. Default 300.
    """
    try:
        import loky
    except ImportError:
        return
    ex = _reusable_executor(max_workers, pool_timeout)
    n = max_workers or loky.cpu_count()
    list(ex.map(int, range(2 * n)))  # force the workers to start


def _map_parallel(fn, arglist, max_workers, parallel, pool_timeout):
    """Run ``fn(*args)`` over ``arglist``, in parallel by default.

    Uses loky's reusable executor if installed, else the stdlib ProcessPoolExecutor.

    Parameters
    ----------
    fn : callable
        The function to map; called as ``fn(*args)``.
    arglist : list of tuple
        Argument tuples, one per point.
    max_workers : int or None
        Process count; None lets the executor choose (~cpu count).
    parallel : bool
        False runs serially in-process.
    pool_timeout : float
        loky idle timeout [s] (see :func:`_reusable_executor`).

    Returns
    -------
    list
        The ``fn`` results, in input order.
    """
    if not parallel:
        return [fn(*args) for args in arglist]
    try:
        ex = _reusable_executor(max_workers, pool_timeout)
        return [f.result() for f in [ex.submit(fn, *args) for args in arglist]]
    except ImportError:
        from concurrent.futures import ProcessPoolExecutor

        with ProcessPoolExecutor(max_workers=max_workers) as ex:
            return [f.result() for f in [ex.submit(fn, *args) for args in arglist]]


def _map_serial_cutoff(fn, arglist, time_cutoff, max_failures):
    """Map ``fn`` over ``arglist`` serially, abandoning the rest once a cutoff is crossed.

    Parameters
    ----------
    fn : callable
        The function to map; called as ``fn(*args)``.
    arglist : list of tuple
        Argument tuples, one per point.
    time_cutoff : float or None
        Abandon once elapsed wall time [s] exceeds this; None disables the time cutoff.
    max_failures : int or None
        Abandon once more than this many points fail; None disables the failure cutoff.
        ``max_failures=0`` bails on the first failure.

    Returns
    -------
    list
        One entry per input point: the ``fn`` result, or ``None`` for points left unrun after the
        cutoff (the caller pads them as skipped-result dicts so the list stays full length).
    """
    t0 = timer()
    out, n_fail, stop = [], 0, False
    for args in arglist:
        if stop:
            out.append(None)
            continue
        r = fn(*args)
        out.append(r)
        if (not r.get("converged", False)) or ("error" in r):
            n_fail += 1
        if (time_cutoff is not None and timer() - t0 > time_cutoff) or \
           (max_failures is not None and n_fail > max_failures):
            stop = True
    return out


def _solve_safe(cfg, scfg, warmstart, keep_weights):
    """Run :func:`solve`, returning an error dict instead of raising (for parallel workers).

    Parameters
    ----------
    cfg : DiffParams
        Problem definition.
    scfg : SolverCfg
        Solver settings.
    warmstart : dict or None
        Optional snapshot to seed the solve.
    keep_weights : bool
        Passed through to :func:`solve`.

    Returns
    -------
    dict
        The :func:`solve` result, or on failure a result dict with ``error`` and ``traceback``
        strings.
    """
    try:
        return solve(cfg, scfg, warmstart=warmstart, keep_weights=keep_weights)
    except Exception as exc:  # noqa: BLE001 (report any failure as data)
        import traceback
        return {"bvalue": 0.0, "converged": False, "n_iter": 0, "n_feval": 0, "X": None,
                "error": f"{type(exc).__name__}: {exc}", "traceback": traceback.format_exc()}


def _cfg_timing_error(cfg: DiffParams):
    """Replicate the C++ ``diff_init*`` timing to reject infeasible geometries before they reach C++.

    An infeasible timing (e.g. ``TE < T_readout``) gives a negative ``N`` or an index outside
    ``[0, N)``, which aborts the process inside Eigen and would take down a parallel batch. Only
    timings that would crash are flagged.

    Parameters
    ----------
    cfg : DiffParams
        Problem definition whose timing is checked.

    Returns
    -------
    str or None
        None if the geometry is safe, else a short reason string.
    """
    dt, TE = cfg.dt, cfg.TE
    if dt <= 0 or TE <= 0:
        return f"non-positive dt/TE (dt={dt}, TE={TE})"
    n_pre = int(np.floor((cfg.T_pre or 0.0) / dt)) if cfg.diff_mode == "preencode" else 0
    N = n_pre + int((TE - cfg.T_readout) / dt) + 1
    if N < cfg.MMT + 3:
        return f"waveform too short: N={N} for MMT={cfg.MMT} (need >= MMT+3; check TE vs T_readout/dt)"
    ind_90_end = n_pre + int(np.ceil(cfg.T_90 / dt))
    if cfg.diff_mode == "conventional":
        ind_180_end = int(np.ceil((TE / 2.0 + cfg.T_180 / 2.0) / dt))
        ind_180_start = ind_90_end + (N - ind_180_end) - 1
    else:
        ind_180_start = n_pre + int(np.floor((TE / 2.0 - cfg.T_180 / 2.0) / dt))
        ind_180_end = n_pre + int(np.ceil((TE / 2.0 + cfg.T_180 / 2.0) / dt))
    for nm, idx in (("ind_90_end", ind_90_end), ("ind_180_start", ind_180_start), ("ind_180_end", ind_180_end)):
        if idx < 0 or idx >= N:
            return (f"{nm}={idx} outside [0,{N}) -- infeasible timing "
                    f"(TE={TE}, T_180={cfg.T_180}, T_readout={cfg.T_readout})")
    return None


def sweep_points(points, *, warmstart: dict = None, keep_weights: bool = True, parallel: bool = True,
                 max_workers: int = None, pool_timeout: float = 300.0, robust: bool = True,
                 time_cutoff: float = None, max_failures: int = None):
    """Solve an explicit list of points in parallel (the list-based version of :func:`sweep`).

    Use this when the points are not a Cartesian grid, e.g. random samples from :func:`sample_points`.
    Points run on the same reusable loky pool as :func:`sweep` / :func:`warm_pool`.

    Parameters
    ----------
    points : list
        Each item is a ``(DiffParams, SolverCfg)`` pair, or a bare ``DiffParams`` (paired with a default
        ``SolverCfg``).
    warmstart, keep_weights, parallel, max_workers, pool_timeout
        As in :func:`sweep`.
    robust : bool, optional
        If True (default), a point that raises returns a result dict with an ``error`` key instead of
        aborting the whole batch. Set False to let exceptions propagate (debugging).
    time_cutoff : float, optional
        Abandon the sweep once its elapsed wall time [s] exceeds this. Serial only (``parallel=False``).
    max_failures : int, optional
        Abandon the sweep once more than this many points have failed (not converged or errored); ``0``
        stops at the first failure. Serial only. Unreached points are returned with
        ``error="abandoned"`` so the list stays full length.

    Returns
    -------
    list of dict
        One result per point, in input order, each also carrying ``"cfg"`` and ``"scfg"``.
    """
    pts = [(p if isinstance(p, tuple) else (p, None)) for p in points]
    pts = [(cfg, scfg or SolverCfg()) for (cfg, scfg) in pts]
    fn = _solve_safe if robust else solve

    early = time_cutoff is not None or max_failures is not None
    if early and parallel:
        msg = "time_cutoff / max_failures require parallel=False (early abandonment is serial)"
        raise ValueError(msg)

    # robust mode: never dispatch timings that would abort the C++ process
    errs = [(_cfg_timing_error(cfg) if robust else None) for (cfg, _) in pts]
    live = [i for i, e in enumerate(errs) if e is None]
    live_args = [(pts[i][0], pts[i][1], warmstart, keep_weights) for i in live]
    if early:
        solved = _map_serial_cutoff(fn, live_args, time_cutoff, max_failures)
    else:
        solved = _map_parallel(fn, live_args, max_workers, parallel, pool_timeout)

    results = [None] * len(pts)
    for i, r in zip(live, solved, strict=True):
        results[i] = r
    for i, e in enumerate(errs):
        if e is not None:
            results[i] = {"bvalue": 0.0, "converged": False, "n_iter": 0, "n_feval": 0, "X": None,
                          "error": f"invalid config: {e}"}
    # pad points skipped by an early cutoff so the list stays full length
    for i in range(len(pts)):
        if results[i] is None:
            results[i] = {"bvalue": 0.0, "converged": False, "n_iter": 0, "n_feval": 0, "X": None,
                          "error": "abandoned"}
    for (cfg, scfg), r in zip(pts, results, strict=True):
        r["cfg"] = cfg
        r["scfg"] = scfg
    return results


def was_abandoned(results) -> bool:
    """Return True if a sweep was cut short by ``time_cutoff`` / ``max_failures``.

    Abandoned points are padded with ``error='abandoned'``::

        res = sweep(cfg, scfg, cfg_grid=grid, parallel=False, time_cutoff=45, max_failures=0)
        if not was_abandoned(res):
            ...  # keep the run

    Parameters
    ----------
    results : list of dict
        The result list returned by :func:`sweep` / :func:`sweep_points`.

    Returns
    -------
    bool
        True if any point carries ``error='abandoned'``.
    """
    return any(r.get("error") == "abandoned" for r in results)


def sample_points(base_cfg: DiffParams, base_scfg: SolverCfg = None, *, cfg_dists=None, scfg_dists=None,
                  n: int = 100, seed: int = None):
    """Draw ``n`` random ``(DiffParams, SolverCfg)`` points to feed :func:`sweep_points`.

    ``cfg_dists`` / ``scfg_dists`` map a field name to a sampler, one of:

    * a callable ``rng -> value``, e.g. log-uniform ``lambda r: 10 ** r.uniform(-3, -1)``;
    * a 2-tuple ``(lo, hi)``: uniform float in ``[lo, hi)``;
    * a list/tuple/array of choices, picked uniformly (use this for ints/strings, e.g. ``MMT``).

    Fields not listed keep the ``base_cfg`` / ``base_scfg`` value; ``rng`` is a
    ``numpy.random.Generator``. A 2-tuple of numbers is a uniform range, so use a 2-element list for
    a two-way choice.

    Parameters
    ----------
    base_cfg, base_scfg : DiffParams, SolverCfg
        The template; unsampled fields are taken from here (``base_scfg`` defaults to ``SolverCfg()``).
    cfg_dists, scfg_dists : dict {field: sampler}, optional
        Per-field samplers as described above.
    n : int, optional
        Number of points to draw. Default 100.
    seed : int, optional
        Seed for the random generator.

    Returns
    -------
    list of tuple
        ``(DiffParams, SolverCfg)`` points.
    """
    rng = np.random.default_rng(seed)
    base_scfg = base_scfg or SolverCfg()

    def draw(dists):
        out = {}
        for field, spec in (dists or {}).items():
            if callable(spec):
                out[field] = spec(rng)
            elif isinstance(spec, tuple) and len(spec) == 2 and all(isinstance(v, (int, float)) for v in spec):
                out[field] = float(rng.uniform(*spec))
            else:  # list/tuple/array of discrete choices
                out[field] = spec[int(rng.integers(len(spec)))]
        return out

    return [(replace(base_cfg, **draw(cfg_dists)), replace(base_scfg, **draw(scfg_dists)))
            for _ in range(n)]


def sweep(base_cfg: DiffParams, base_scfg: SolverCfg = None, *, cfg_grid=None, scfg_grid=None,
          warmstart: dict = None, keep_weights: bool = True, parallel: bool = True, max_workers: int = None,
          pool_timeout: float = 300.0, time_cutoff: float = None, max_failures: int = None):
    """Solve the cross product of ``cfg_grid`` x ``scfg_grid`` in parallel, one process per point.

    ``DiffParams`` and ``SolverCfg`` field names are disjoint, so one call can vary either or both::

        sweep(cfg, scfg, cfg_grid={"TE": np.linspace(50e-3, 100e-3, 8)})
        sweep(cfg, scfg, scfg_grid={"ils_tol": [0.05, 0.1, 0.2]})
        sweep(cfg, scfg, cfg_grid={"TE": tes}, scfg_grid={"gamma_x": [1.4, 1.6]})

    Parameters
    ----------
    base_cfg, base_scfg : DiffParams, SolverCfg
        The fixed problem/solver; ``base_scfg`` defaults to ``SolverCfg()``.
    cfg_grid, scfg_grid : dict {field: [values]}, optional
        Field overrides applied via ``replace``; ``None`` => that object is not varied.
    warmstart : dict, optional
        One snapshot fed to every point (resized per grid), e.g. from a representative solve. Most
        useful when the points are close to each other.
    keep_weights : bool, optional
        As in :func:`solve`.
    parallel : bool, optional
        False runs serially in-process (debugging / tiny sweeps).
    max_workers : int, optional
        Process count; None lets the executor choose (~cpu count). Uses loky's reusable pool if
        installed, else the stdlib ProcessPoolExecutor.
    pool_timeout : float, optional
        loky only: idle seconds before the worker pool is reaped (default 300; loky's own 10 s default
        reaps it between notebook cells). Use the same value across calls and with :func:`warm_pool`.
    time_cutoff, max_failures : optional
        Early abandonment, as in :func:`sweep_points` (serial only).

    Returns
    -------
    list of dict
        One result per point (cfg_grid outer, scfg_grid inner), each also carrying ``"cfg"`` and
        ``"scfg"``.
    """
    base_scfg = base_scfg or SolverCfg()
    points = _grid_points(base_cfg, base_scfg, cfg_grid, scfg_grid)
    return sweep_points(points, warmstart=warmstart, keep_weights=keep_weights, parallel=parallel,
                        max_workers=max_workers, pool_timeout=pool_timeout,
                        time_cutoff=time_cutoff, max_failures=max_failures)

