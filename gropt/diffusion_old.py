from dataclasses import dataclass, replace
from typing import Optional

import gropt
import numpy as np


def diff_min_TE(params, target_bval=0, TE0=10e-3, TE1=120e-3, stop_dt=0.1e-3, **kwargs):

    if target_bval == 0:
        target_bval = params['bvalue']

    stop_dt = max(stop_dt, 2 * params['dt'])

    min_TE = params['T_90'] + params['T_180'] + params['T_readout']
    min_TE = max(min_TE, 2 * (params['T_180'] + params['T_readout']))

    TE0 += min_TE

    result0 = diff_solve_TE(TE0, params, bval_min=target_bval / 2, **kwargs)
    result1 = diff_solve_TE(TE1, params, bval_min=target_bval / 2, **kwargs)

    for _i in range(3):
        if not result1.converged or result1.bvalue < target_bval:
            result1 = diff_solve_TE(TE1, params, bval_min=4 ** (_i + 1) * target_bval, **kwargs)

    if not result1.converged or result1.bvalue < target_bval:
        print(
            f'ERROR: Max TE {TE1} did not converge or was too low bvalue: {result1.converged}  {result1.bvalue}'
        )

    print('Searching TE: ', end='', flush=True)
    for _i in range(10):
        _TE = TE0 + (TE1 - TE0) / 2
        _result = diff_solve_TE(_TE, params, bval_min=target_bval / 2, **kwargs)

        if not _result.converged:
            print(f'!{_TE * 1000:.2f}', end=' ', flush=True)
        else:
            print(f'{_TE * 1000:.2f}', end=' ', flush=True)

        if _result.converged:
            if _result.bvalue > target_bval:
                TE1 = _TE
                result1 = _result
            else:
                TE0 = _TE
                result0 = _result
        else:
            TE0 = _TE
            result0 = _result

    print('Done!', flush=True)

    if result0.converged and result0.bvalue > target_bval:
        return TE0, result0
    if _result.converged and _result.bvalue > target_bval:
        return _TE, _result
    if result1.converged and result1.bvalue > target_bval:
        return TE1, result1

    return None, None


def diff_solve_TE(TE, params, bval_min=100.0, refine=False, **kwargs):
    result = _diff_solve_TE(TE, params, bval_min=bval_min, **kwargs)

    if refine:
        if not result.converged:
            result2 = _diff_solve_TE(TE, params, bval_min=bval_min / 2, **kwargs)
        else:
            bval0 = result.bvalue
            if bval0 < bval_min:
                result2 = _diff_solve_TE(TE, params, bval_min=0.8 * bval0, **kwargs)
            else:
                result2 = _diff_solve_TE(TE, params, bval_min=0.8 * bval0, **kwargs)

        if (result2.converged and not result.converged) or (
            result2.converged and result2.bvalue > result.bvalue
        ):
            result = result2

    return result


def _diff_solve_TE(
    TE,
    params,
    bval_min=100.0,
    dt=0.0,
    obj_patience=30,
    ils_max_iter=24,
    moment_tol=1e-5,
    bval_scale=1.02,
    **kwargs,
):

    if dt == 0:
        dt = params['dt']

    if 'diff_mode' in params:
        diff_mode = params['diff_mode']
    else:
        diff_mode = 'gropt'

    # Set up the GrOpt problem
    gparams = gropt.GroptParams()

    start_idx = 0
    if diff_mode == 'gropt':
        gparams.diff_init(
            dt=dt,
            TE=TE,
            T_90=params['T_90'],
            T_180=params['T_180'],
            T_readout=params['T_readout'],
        )
    elif diff_mode == 'conventional':
        gparams.diff_init_deadtime(
            dt=dt,
            TE=TE,
            T_90=params['T_90'],
            T_180=params['T_180'],
            T_readout=params['T_readout'],
        )
    elif diff_mode == 'preencode':
        start_idx = gparams.diff_init_preencode(
            dt=dt,
            TE=TE,
            T_90=params['T_90'],
            T_180=params['T_180'],
            T_readout=params['T_readout'],
            T_pre=params['T_pre'],
        )
        params['start_idx'] = start_idx
    else:
        msg = f'Unknown diff_mode: {diff_mode}'
        raise ValueError(msg)

    gparams.add_gmax(params['gmax'])
    gparams.add_smax(params['smax'])
    for _i_moment in range(params['MMT'] + 1):
        gparams.add_moment(_i_moment, 0.0, start_idx=start_idx, tol=moment_tol)

    if 'pns_lim' in params:
        if 'pns_params' in params:
            gparams.add_SAFE(params['pns_lim'], safe_params=params['pns_params'])
        else:
            pns_params, cns_params = gropt.get_random_safe_params()
            params['pns_params'] = pns_params
            gparams.add_SAFE(params['pns_lim'], safe_params=pns_params)

    if 'cns_lim' in params:
        if 'cns_params' in params:
            gparams.add_SAFE(params['cns_lim'], safe_params=params['cns_params'])
        else:
            pns_params, cns_params = gropt.get_random_safe_params()
            params['cns_params'] = cns_params
            gparams.add_SAFE(params['cns_lim'], safe_params=cns_params)

    if 'concomitant' in params:
        gparams.add_concomitant(start_idx=start_idx, weight_mod=1e4)

    if 'eddy_lam' in params:
        gparams.add_eddy(params['eddy_lam'])

    gparams.add_bvalue(bval_min, mode='minval_max', start_idx0=start_idx, max_scale=bval_scale)

    gparams.prepare()

    result = diff_solve(gparams, obj_patience=obj_patience, ils_max_iter=ils_max_iter)

    return result


def diff_solve(gparams, obj_patience=30, ils_max_iter=30):
    solver = gropt.SolverGroptSDMM()
    solver.set_general_params(max_feval=200000, max_iter=20000, gamma_x=1.6, obj_patience=obj_patience)
    solver.set_ils_params(ils_max_iter=ils_max_iter, ils_tol=1e-12, ils_sigma=0.0001, ils_tik_lam=0.0001)
    solver.set_sdmm_params(rw_interval=16, grw_interval=41)
    result = solver.solve(gparams)
    return result


def diff_build_gparams(
    TE,
    params,
    bval_min=100.0,
    dt=0.0,
    moment_tol=1e-5,
    bval_scale=1.02,
    bval_mode = 'minval_max',
):

    if dt == 0:
        dt = params['dt']

    if 'diff_mode' in params:
        diff_mode = params['diff_mode']
    else:
        diff_mode = 'gropt'

    # Set up the GrOpt problem
    gparams = gropt.GroptParams()

    start_idx = 0
    if diff_mode == 'gropt':
        gparams.diff_init(
            dt=dt,
            TE=TE,
            T_90=params['T_90'],
            T_180=params['T_180'],
            T_readout=params['T_readout'],
        )
    elif diff_mode == 'conventional':
        gparams.diff_init_deadtime(
            dt=dt,
            TE=TE,
            T_90=params['T_90'],
            T_180=params['T_180'],
            T_readout=params['T_readout'],
        )
    elif diff_mode == 'preencode':
        start_idx = gparams.diff_init_preencode(
            dt=dt,
            TE=TE,
            T_90=params['T_90'],
            T_180=params['T_180'],
            T_readout=params['T_readout'],
            T_pre=params['T_pre'],
        )
        params['start_idx'] = start_idx
    else:
        msg = f'Unknown diff_mode: {diff_mode}'
        raise ValueError(msg)

    gparams.add_gmax(params['gmax'])
    gparams.add_smax(params['smax'])
    for _i_moment in range(params['MMT'] + 1):
        gparams.add_moment(_i_moment, 0.0, start_idx=start_idx, tol=moment_tol)

    if 'pns_lim' in params:
        if 'pns_params' in params:
            gparams.add_SAFE(params['pns_lim'], safe_params=params['pns_params'])
        else:
            pns_params, cns_params = gropt.get_random_safe_params()
            params['pns_params'] = pns_params
            gparams.add_SAFE(params['pns_lim'], safe_params=pns_params)

    if 'cns_lim' in params:
        if 'cns_params' in params:
            gparams.add_SAFE(params['cns_lim'], safe_params=params['cns_params'])
        else:
            pns_params, cns_params = gropt.get_random_safe_params()
            params['cns_params'] = cns_params
            gparams.add_SAFE(params['cns_lim'], safe_params=cns_params)

    if 'concomitant' in params:
        gparams.add_concomitant(start_idx=start_idx, weight_mod=1e4)

    if 'eddy_lam' in params:
        gparams.add_eddy(params['eddy_lam'])

    if bval_mode == 'obj':
        gparams.add_bvalue(as_objective=True, start_idx0=start_idx, weight_mod=bval_scale)
    else:
        gparams.add_bvalue(bval_min, mode=bval_mode, start_idx0=start_idx, max_scale=bval_scale)

    gparams.prepare()

    return gparams


# ===========================================================================
# Dataclass-based config 
# ===========================================================================
@dataclass(frozen=True)
class DiffParams:
    # Timing
    TE: float
    T_90: float
    T_180: float
    T_readout: float
    dt: float = 400e-6
    # hardware
    gmax: float = 0.08  # [T/m]
    smax: float = 200.0  # [T/m/s]
    # diffusion
    MMT: int = 0  # null M0..M_MMT
    bvalue: float = 1000.0  # target [s/mm^2]
    # mode / structure
    diff_mode: str = "gropt"  # 'gropt' | 'conventional' | 'preencode'
    T_pre: Optional[float] = None  # preencode only
    # optional constraints (None / False => off)
    concomitant: bool = False
    eddy_lam: Optional[float] = None
    pns_lim: Optional[float] = None
    cns_lim: Optional[float] = None
    safe_seed: int = 42                 # regenerate SAFE dicts in-worker (hashable + reproducible)
    # problem-shaping knobs
    moment_tol: float = 1e-5
    moment_project: bool = True
    bval_mode: str = "obj"              # 'obj' (maximize) or a constraint mode e.g. 'minval_max'
    bval_min: float = 100.0             # only used when bval_mode != 'obj'
    bval_scale: float = 1.02

    def structural_key(self):
        """Discrete signature for bucketing warm-start matches (exact match required)."""
        return (self.diff_mode, self.MMT, self.concomitant, self.eddy_lam is not None,
                self.pns_lim is not None, self.cns_lim is not None, self.bval_mode)


def build_gparams(cfg: DiffParams):
    """Build a GroptParams from a DiffParams. Returns (gparams, start_idx); does NOT mutate cfg."""
    gp = gropt.GroptParams()

    start_idx = 0
    if cfg.diff_mode == "gropt":
        gp.diff_init(dt=cfg.dt, TE=cfg.TE, T_90=cfg.T_90, T_180=cfg.T_180, T_readout=cfg.T_readout)
    elif cfg.diff_mode == "conventional":
        gp.diff_init_deadtime(dt=cfg.dt, TE=cfg.TE, T_90=cfg.T_90, T_180=cfg.T_180, T_readout=cfg.T_readout)
    elif cfg.diff_mode == "preencode":
        start_idx = gp.diff_init_preencode(dt=cfg.dt, TE=cfg.TE, T_90=cfg.T_90, T_180=cfg.T_180,
                                           T_readout=cfg.T_readout, T_pre=cfg.T_pre)
    else:
        raise ValueError(f"Unknown diff_mode: {cfg.diff_mode}")

    gp.add_gmax(cfg.gmax)
    gp.add_smax(cfg.smax)
    for m in range(cfg.MMT + 1):
        gp.add_moment(m, 0.0, start_idx=start_idx, tol=cfg.moment_tol, project=cfg.moment_project)

    if cfg.pns_lim is not None:
        pns, _ = gropt.get_random_safe_params(cfg.safe_seed)
        gp.add_SAFE(cfg.pns_lim, safe_params=pns)
    if cfg.cns_lim is not None:
        _, cns = gropt.get_random_safe_params(cfg.safe_seed)
        gp.add_SAFE(cfg.cns_lim, safe_params=cns)
    if cfg.concomitant:
        gp.add_concomitant(start_idx=start_idx, project=True)
    if cfg.eddy_lam is not None:
        gp.add_eddy(cfg.eddy_lam)

    if cfg.bval_mode == "obj":
        gp.add_bvalue(as_objective=True, start_idx0=start_idx, weight_mod=cfg.bval_scale)
        gp.normalize_obj = True
    else:
        gp.add_bvalue(cfg.bval_min, mode=cfg.bval_mode, start_idx0=start_idx, max_scale=cfg.bval_scale)

    gp.prepare()
    return gp, start_idx


def add_basin_anchor(gparams, cfg, window_time=1e-3, eps_factor=0.07, same_sign=False):
    """Add the diffusion basin-orientation constraint (GroptParams.add_diff_basin), pulling gmax
    from cfg so you only pass the window/eps tuning. cfg may be a DiffParams or a params dict.

    Keep eps_factor small (it's inactive at the optimum) for single-solve basin control, and keep
    its sign convention consistent with any X0 seed orientation. same_sign=False (default) forces the
    opposite-sign flip basin; same_sign=True forces the no-flip (local) basin.
    """
    gmax = cfg.gmax if hasattr(cfg, "gmax") else cfg["gmax"]
    gparams.add_diff_basin(window_time, eps_factor, gmax, same_sign=same_sign)


def solve_one(cfg: DiffParams, solver_kwargs=None):
    """loky-friendly worker: frozen cfg in, plain picklable data out (build gparams here).

    solver_kwargs : dict | None
        Arbitrary Solver attribute overrides applied via setattr, e.g.
        {"max_iter": 4000, "ils_tol": 0.01, "gamma_x": 1.4, "obj_patience": 10}.
        Kept OUT of DiffParams so the config stays hashable/picklable. Only the
        read/write Solver attributes work this way; method-only knobs (e.g.
        set_sdmm_params: rw_interval, grw_interval) still need their own call.
    """
    gp, _ = build_gparams(cfg)
    solver = gropt.SolverGroptSDMM()
    # placeholder defaults; override/extend via solver_kwargs
    solver.max_iter = 4000
    solver.obj_patience = 10
    solver.ils_max_iter = 20
    for k, v in (solver_kwargs or {}).items():
        if not hasattr(solver, k):  # guard: setattr has no validation, so a typo would silently no-op
            raise AttributeError(f"SolverGroptSDMM has no attribute '{k}' -- check solver_kwargs for typos")
        setattr(solver, k, v)
    r = solver.solve(gp)
    return dict(TE=cfg.TE, bvalue=float(r.bvalue), converged=bool(r.converged), X=np.asarray(r.X))

