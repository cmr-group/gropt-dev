import matplotlib.pyplot as plt
import numpy as np


def get_moments(g, dt, inv_vec=None, start_idx=0, scale_to_one=True):
    if g.squeeze().ndim == 2:
        g = g[0]  # TODO: 3-axis case, right now just assumes 1 axis

    Nm = 5
    tt = np.arange(g.size) * dt
    _m = np.zeros((Nm, g.size))

    for mm in range(Nm):
        _m[mm] = tt**mm

    if inv_vec is not None:
        mm = dt * _m * (g * inv_vec)[np.newaxis, :]
    else:
        mm = dt * _m * g[np.newaxis, :]

    mm[:, :start_idx] = 0

    moments = [np.cumsum(x) for x in mm]

    if scale_to_one:
        moments = [x / np.abs(x).max() for x in moments]

    return moments


def get_concomitant(g, dt, inv_vec, start_idx=0):
    g_start = g[start_idx:]
    inv_vec_start = inv_vec[start_idx:]
    pos = np.sum(dt * g_start[inv_vec_start > 0] ** 2.0)
    neg = np.sum(dt * g_start[inv_vec_start < 0] ** 2.0)

    ratio = pos / neg
    if ratio < 1:
        ratio = 1 / ratio

    return ratio


def get_bval(g, dt, inv_vec=None, TE=0, start_idx=0):
    if g.squeeze().ndim == 2:
        g = g[0]  # TODO: 3-axis case, right now just assumes 1 axis

    if inv_vec is None:
        inv_vec = np.ones(g.size)
        tINV = int(np.floor(TE / dt / 2.0))
        inv_vec[tINV:] = -1

    GAMMA = 42.58e3

    Gt = 0
    bval = 0
    for i in range(start_idx, g.size):
        Gt += inv_vec[i] * g[i] * dt
        bval += Gt * Gt * dt

    bval *= (GAMMA * 2 * np.pi) ** 2

    return bval


def plot_diff(*args, mode='diff', **kwargs):

    plot_waves(*args, mode=mode, **kwargs)


def plot_waves(
    g,
    dt,
    inv_vec=None,
    TE=0,
    gmax=0,
    smax=0,
    start_idx=0,
    eddy_lam=0.0,
    plot_eddy=False,
    safe_params=None,
    MMT=4,
    stim=None,
    stim_vec=None,
    stim_lim=0,
    N_cols=2,
    mode='regular',
    params={},
):

    if start_idx == 0:
        start_idx = params.get('start_idx', 0)
    if eddy_lam == 0:
        eddy_lam = params.get('eddy_lam', 0.0)
    if stim_vec is None:
        stim_vec = params.get('stim_vec', None)
    if stim_lim == 0:
        stim_lim = params.get('stim_lim', 0.0)

    if inv_vec is None:
        inv_vec = np.ones(g.size)
        tINV = start_idx + int(np.floor(TE / dt / 2.0))
        inv_vec[tINV:] = -1

    tt_ms = np.arange(g.size) * dt * 1e3

    to_plot = ['gradient', 'slew', 'moments']
    if stim is not None:
        to_plot.append('stim')
    if plot_eddy or eddy_lam > 0:
        to_plot.append('eddy')

    N_plots = len(to_plot)
    N_rows = (N_plots + N_cols - 1) // N_cols
    f, axarr = plt.subplots(N_rows, N_cols, squeeze=False, figsize=(10, N_rows * 3.0), layout='tight')

    # Diffusion Title String
    # ======================================
    if mode == 'diff':
        label = ''

        if TE > 0:
            label += f'TE: {1000 * TE:.2f} ms  ---  '

        bval = get_bval(g, dt, inv_vec, start_idx=start_idx)
        label += f'b-value: {bval:.0f} $mm^2/s$  ---  '

        c_ratio = get_concomitant(g, dt, inv_vec, start_idx=start_idx)
        label += f'concomitant ratio: {c_ratio:.2f}'

        f.suptitle(label)

    for i_ax, ax in enumerate(axarr.flatten()):
        plot_type = to_plot[i_ax] if i_ax < N_plots else None
        if plot_type is None:
            ax.set_visible(False)
        elif plot_type == 'gradient':
            ax.axhline(linestyle='--', color='0.7')
            if gmax == 0:
                gmax = params.get('gmax', 0)
            if gmax > 0:
                ax.axhline(1000 * gmax, linestyle=':', color='r', alpha=0.5)
                ax.axhline(-1000 * gmax, linestyle=':', color='r', alpha=0.5)
            ax.plot(tt_ms, g * 1000)
            ax.set_title('Gradient')
            ax.set_xlabel('t [ms]')
            ax.set_ylabel('G [mT/m]')
        elif plot_type == 'slew':
            ax.axhline(linestyle='--', color='0.7')
            if smax == 0:
                smax = params.get('smax', 0)
            if smax > 0:
                ax.axhline(smax, linestyle=':', color='r', alpha=0.5)
                ax.axhline(-smax, linestyle=':', color='r', alpha=0.5)
            ax.plot(tt_ms[:-1], np.diff(g) / dt)
            ax.set_title('Slew')
            ax.set_xlabel('t [ms]')
            ax.set_ylabel('dG/dt [T/m/s]')
        elif plot_type == 'moments':
            mm = get_moments(g, dt, inv_vec, start_idx=start_idx)

            ax.axhline(linestyle='--', color='0.7')
            for im in range(params.get('MMT', 4) + 1):
                ax.plot(tt_ms, mm[im], label=f'{im}')
            ax.legend(loc='upper left')
            ax.set_title('Moments')
            ax.set_xlabel('t [ms]')
            ax.set_ylabel('Moment [a.u.]')
        elif plot_type == 'stim':
            ax.plot(tt_ms, stim)
            if stim_vec is not None:
                ax.plot(tt_ms, stim_vec, linestyle=':', color='r', alpha=0.5)
            if stim_lim > 0:
                ax.axhline(stim_lim, linestyle=':', color='r', alpha=0.5)
            ax.axhline(linestyle='--', color='0.7')

            ax.set_title('SAFE PNS')
            ax.set_xlabel('t [ms]')
        elif plot_type == 'eddy':
            all_lam = np.linspace(0.1, 120, 1000)
            all_e = []
            for lam in all_lam:
                _lam = lam * 1.0e-3
                r = np.diff(np.exp(-np.arange(g.size + 1) * dt / _lam))[::-1]
                all_e.append(100 * r @ g)
            all_e = np.array(all_e)

            if all_e.min() < -all_e.max():
                all_e = -all_e

            ax.axhline(linestyle='--', color='0.7')
            ax.plot(all_lam, all_e)

            if eddy_lam > 0:
                ax.axvline(eddy_lam * 1e3, linestyle='--', color='r', alpha=0.7)

            ax.set_title('Eddy Spectrum')
            ax.set_xlabel('lambda [ms]')
            ax.set_ylabel('eddy spectrum [a.u.]')

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    plt.show()
