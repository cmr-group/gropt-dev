import gropt
import numpy as np


def diff_min_TE(params, target_bval=0, TE0=10e-3, TE1=120e-3, stop_dt=0.1e-3):

    if target_bval == 0:
        target_bval = params['bvalue']

    stop_dt = max(stop_dt, 2 * params['dt'])

    min_TE = params['T_90'] + params['T_180'] + params['T_readout']
    min_TE = max(min_TE, 2 * (params['T_180'] + params['T_readout']))

    TE0 += min_TE

    result0 = diff_solve_TE(TE0, params, bval_min=target_bval / 2)
    result1 = diff_solve_TE(TE1, params, bval_min=target_bval / 2)

    if not result1.converged or result1.bvalue < target_bval:
        result1 = diff_solve_TE(TE1, params, bval_min=4 * target_bval)
    if not result1.converged or result1.bvalue < target_bval:
        print(
            f'ERROR: Max TE {TE1} did not converge or was too low bvalue: {result1.converged}  {result1.bvalue}'
        )

    print('Searching TE: ', end='', flush=True)
    for _i in range(10):
        _TE = TE0 + (TE1 - TE0) / 2
        _result = diff_solve_TE(_TE, params, bval_min=target_bval / 2)

        # print(
        #     f'TE: {_TE * 1000:.2f}   [{TE0 * 1000:.2f} {TE1 * 1000:.2f}]  {_result.converged}  bval: {_result.bvalue:.2f}'
        # )
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


def _diff_solve_TE(TE, params, bval_min=100.0, dt=0.0, extra_iters=1000, ils_max_iter=30):

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
        gparams.add_moment(_i_moment, 0.0, start_idx=start_idx)

    if 'stim_lim' in params:
        gparams.add_SAFE(params['stim_lim'])

    if 'concomitant' in params:
        gparams.add_concomitant(start_idx=start_idx)

    if 'eddy_lam' in params:
        gparams.add_eddy(params['eddy_lam'])

    gparams.add_bvalue(bval_min, mode='minval_max', start_idx0=start_idx)

    gparams.prepare()

    result = diff_solve(gparams, extra_iters=extra_iters, ils_max_iter=ils_max_iter)

    return result


def diff_solve(gparams, extra_iters=2000, ils_max_iter=30):
    solver = gropt.SolverGroptSDMM()
    solver.set_general_params(max_feval=200000, max_iter=20000, gamma_x=1.5, extra_iters=extra_iters)
    solver.set_ils_params(ils_max_iter=ils_max_iter, ils_tol=1e-12, ils_sigma=0.0001, ils_tik_lam=0.0001)
    result = solver.solve(gparams)
    return result
