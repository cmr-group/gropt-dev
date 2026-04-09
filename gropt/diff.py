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
    extra_iters=1000,
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

    if 'tv_lam' in params:
        gparams.add_TV(params['tv_lam'], weight_mod=params['tv_weight'] if 'tv_weight' in params else 5)

    if 'identity_lam' in params:
        gparams.add_obj_identity(params['identity_lam'])
    
    if 'acoustic' in params:
        freqs = params['acoustic']['freqs']
        bws = params['acoustic'].get('bws', np.zeros_like(freqs))
        if len(freqs) > 0:
            gparams.add_acoustic(freqs, bws, weight_mod=0.5, bw_scale=0.6)

    gparams.add_bvalue(bval_min, mode='minval_max', start_idx0=start_idx, max_scale=bval_scale)

    gparams.prepare()

    result = diff_solve(gparams, extra_iters=extra_iters, ils_max_iter=ils_max_iter)

    return result


def diff_solve(gparams, extra_iters=2000, ils_max_iter=30):
    solver = gropt.SolverGroptSDMM()
    solver.set_general_params(max_feval=200000, max_iter=20000, gamma_x=1.6, extra_iters=extra_iters)
    solver.set_ils_params(ils_max_iter=ils_max_iter, ils_tol=1e-12, ils_sigma=0.0001, ils_tik_lam=0.0001)
    solver.set_sdmm_params(rw_interval=16, grw_interval=41)
    result = solver.solve(gparams)
    return result
