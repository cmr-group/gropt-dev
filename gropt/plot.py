import matplotlib.pyplot as plt
import numpy as np


def get_moments(g, dt, inv_vec=None, scale_to_one=True):
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

    moments = [np.cumsum(x) for x in mm]

    if scale_to_one:
        moments = [x / np.abs(x).max() for x in moments]

    return moments


def get_concomitant(g, dt, inv_vec):
    pos = np.sum(dt * g[inv_vec > 0] ** 2.0)
    neg = np.sum(dt * g[inv_vec < 0] ** 2.0)

    ratio = pos / neg
    if ratio < 1:
        ratio = 1 / ratio

    return ratio


def get_bval(g, dt, inv_vec=None, TE=0):
    if g.squeeze().ndim == 2:
        g = g[0]  # TODO: 3-axis case, right now just assumes 1 axis

    if inv_vec is None:
        inv_vec = np.ones(g.size)
        tINV = int(np.floor(TE / dt / 2.0))
        inv_vec[tINV:] = -1

    GAMMA = 42.58e3

    Gt = 0
    bval = 0
    for i in range(g.size):
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
    safe_params=None,
    MMT=4,
    stim=None,
    stim_vec=None,
    stim_lim=0,
    mode='regular',
):

    if inv_vec is None:
        inv_vec = np.ones(g.size)
        tINV = int(np.floor(TE / dt / 2.0))
        inv_vec[tINV:] = -1

    tt_ms = np.arange(g.size) * dt * 1e3

    N_rows = 2
    N_cols = 2
    f, axarr = plt.subplots(N_rows, N_cols, squeeze=False, figsize=(10, N_rows * 3.0), layout='constrained')

    i_row = 0
    i_col = 0

    if mode == 'diff':
        label = ''

        if TE > 0:
            label += f'TE: {1000 * TE} ms  ---  '

        bval = get_bval(g, dt, inv_vec)
        label += f'b-value: {bval:.0f} $mm^2/s$  ---  '

        c_ratio = get_concomitant(g, dt, inv_vec)
        label += f'concomitant ratio: {c_ratio:.1f}'

        f.suptitle(label)

    axarr[i_row, i_col].axhline(linestyle='--', color='0.7')
    if gmax > 0:
        axarr[i_row, i_col].axhline(1000 * gmax, linestyle=':', color='r', alpha=0.5)
        axarr[i_row, i_col].axhline(-1000 * gmax, linestyle=':', color='r', alpha=0.5)
    axarr[i_row, i_col].plot(tt_ms, g * 1000)
    axarr[i_row, i_col].set_title('Gradient')
    axarr[i_row, i_col].set_xlabel('t [ms]')

    i_row = 0
    i_col = 1

    axarr[i_row, i_col].axhline(linestyle='--', color='0.7')
    if smax > 0:
        axarr[i_row, i_col].axhline(smax, linestyle=':', color='r', alpha=0.5)
        axarr[i_row, i_col].axhline(-smax, linestyle=':', color='r', alpha=0.5)
    axarr[i_row, i_col].plot(tt_ms[:-1], np.diff(g) / dt)
    axarr[i_row, i_col].set_title('Slew')
    axarr[i_row, i_col].set_xlabel('t [ms]')

    i_row = 1
    i_col = 0

    mm = get_moments(g, dt, inv_vec)

    axarr[i_row, i_col].axhline(linestyle='--', color='0.7')
    for im in range(5):
        axarr[i_row, i_col].plot(tt_ms, mm[im], label=f'{im}')
    axarr[i_row, i_col].legend(loc='upper left')
    axarr[i_row, i_col].set_title('Moments')
    axarr[i_row, i_col].set_xlabel('t [ms]')

    i_row = 1
    i_col = 1

    if stim is not None:
        axarr[i_row, i_col].plot(tt_ms, stim)

        if stim_vec is not None:
            axarr[i_row, i_col].plot(tt_ms, stim_vec, linestyle=':', color='r', alpha=0.5)
        if stim_lim > 0:
            axarr[i_row, i_col].axhline(stim_lim, linestyle=':', color='r', alpha=0.5)

        axarr[i_row, i_col].axhline(linestyle='--', color='0.7')

    axarr[i_row, i_col].set_title('SAFE PNS')
    axarr[i_row, i_col].set_xlabel('t [ms]')

    plt.show()
