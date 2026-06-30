"""Interactive / diagnostic plots for the gropt solver's debug history.

Everything here consumes the dict returned by ``solver.get_debug()`` (populated
only when ``solver.extra_debug = True`` before ``solve()``). Plotting backends
(matplotlib / plotly / ipywidgets) are imported lazily inside each function so
that ``import gropt`` stays light and these heavy deps stay optional.

Typical use in a notebook::

    solver.extra_debug = True
    result = gropt.solve(gparams, solver)
    dbg = solver.get_debug()

    from gropt import debug_plotting as dbp
    dbp.scrub_hist_X(dbg, Naxis=gparams.Naxis, dt=gparams.dt)  # plotly scrubber
    dbp.plot_convergence(dbg)                                  # solve dashboard
"""

import numpy as np


# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #
def _as2d(seq):
    """Per-iteration / per-operator history -> (n_iter, n_op) float array, or None."""
    if seq is None or len(seq) == 0:
        return None
    try:
        a = np.array(seq, dtype=float)
    except Exception:
        return None
    if a.size == 0:
        return None
    if a.ndim == 1:
        a = a[:, None]
    return a


def _as1d(seq):
    """1-D history -> float array, or None if empty."""
    if seq is None or len(seq) == 0:
        return None
    a = np.asarray(seq, dtype=float)
    return a if a.size else None


def _op_label(i, op_names):
    if op_names is not None and i < len(op_names):
        return op_names[i]
    return f"op{i}"


def _reshape_hist_X(debug, Naxis):
    """hist_X (list of flat vectors) -> (n_iter, Naxis, N) array. Raises if empty."""
    hX = debug.get("hist_X", [])
    if hX is None or len(hX) == 0:
        raise ValueError("hist_X is empty -- set solver.extra_debug = True before solve()")
    n = len(hX)
    N = len(hX[0]) // Naxis
    return np.asarray(hX, dtype=float).reshape(n, Naxis, N), n, N


# --------------------------------------------------------------------------- #
# waveform-history scrubbers
# --------------------------------------------------------------------------- #
def scrub_hist_X(debug, Naxis=1, dt=None, show_slew=True, subsample=1):
    """Interactive plotly slider over the gradient-waveform history (hist_X).

    All frames are embedded client-side, so the slider is pure JS -- no kernel
    round-trip per frame. This is the one to use in VS Code notebooks: it won't
    lag, freeze, or disturb your matplotlib figures.

    Parameters
    ----------
    debug : dict
        Output of ``solver.get_debug()``.
    Naxis : int
        Number of gradient axes (use ``gparams.Naxis``).
    dt : float or None
        Raster time [s]. If given, the x-axis is time [ms] and a slew panel is
        shown; otherwise the x-axis is sample index and slew is skipped.
    show_slew : bool
        Show the slew (first-difference) panel below the gradient panel.
    subsample : int
        Keep every ``subsample``-th iteration. Increase to shrink the notebook
        output for very long runs (all frames are serialized into the cell).

    Returns
    -------
    plotly.graph_objects.Figure
        The figure. Call this as the last expression in a cell so the notebook
        auto-displays it once (it does NOT call fig.show(), which double-renders
        in VS Code).

    Notes
    -----
    Requires ``plotly`` (and ``nbformat`` for the VS Code renderer)::

        pixi add --feature dev plotly nbformat
    """
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except ImportError as e:
        raise ImportError(
            "scrub_hist_X needs plotly. Install with `pixi add --feature dev plotly nbformat`."
        ) from e

    g, n, N = _reshape_hist_X(debug, Naxis)
    if subsample > 1:
        g = g[::subsample]
    idx = np.arange(0, n, subsample)  # original iteration numbers of the kept frames
    nf = g.shape[0]

    t = np.arange(N) * dt * 1e3 if dt else np.arange(N).astype(float)
    xl = "time [ms]" if dt else "sample"
    slew = show_slew and (dt is not None) and N > 1
    if slew:
        s = np.diff(g, axis=2) / dt
        ts = t[:-1]

    bval = _as1d(debug.get("hist_bvalue"))
    feas = _as1d(debug.get("hist_all_feas"))
    best = int(debug.get("best_feasible_iter", -1))

    gpad = 0.05 * (g.max() - g.min() + 1e-12)
    grange = [g.min() - gpad, g.max() + gpad]
    if slew:
        spad = 0.05 * (s.max() - s.min() + 1e-12)
        srange = [s.min() - spad, s.max() + spad]

    def title(k):
        it = int(idx[k])
        p = [f"iter {it}/{n - 1}"]
        if bval is not None and it < len(bval):
            p.append(f"b={bval[it]:.1f}")
        if feas is not None and it < len(feas):
            p.append("feasible" if feas[it] else "INFEASIBLE")
        if it == best:
            p.append("&#8592; returned")
        return "   ".join(p)

    rows = 2 if slew else 1
    fig = make_subplots(rows=rows, cols=1, shared_xaxes=True, vertical_spacing=0.08)

    # frame-0 traces define the trace order that frames map onto by index
    for j in range(Naxis):
        fig.add_trace(go.Scatter(x=t, y=g[0, j], mode="lines", name=f"axis {j}"), row=1, col=1)
    if slew:
        for j in range(Naxis):
            fig.add_trace(go.Scatter(x=ts, y=s[0, j], mode="lines", showlegend=False), row=2, col=1)

    frames = []
    for k in range(nf):
        data = [go.Scatter(x=t, y=g[k, j]) for j in range(Naxis)]
        if slew:
            data += [go.Scatter(x=ts, y=s[k, j]) for j in range(Naxis)]
        frames.append(go.Frame(data=data, name=str(k), layout=go.Layout(title_text=title(k))))
    fig.frames = frames

    steps = [
        dict(
            method="animate",
            label=str(int(idx[k])),
            args=[[str(k)], dict(mode="immediate", frame=dict(duration=0, redraw=True),
                                 transition=dict(duration=0))],
        )
        for k in range(nf)
    ]
    fig.update_layout(
        height=300 * rows,
        title_text=title(0),
        margin=dict(t=60),
        sliders=[dict(active=0, steps=steps, currentvalue=dict(prefix="iter "))],
    )
    fig.update_yaxes(title_text="gradient", range=grange, row=1, col=1)
    if slew:
        fig.update_yaxes(title_text="slew", range=srange, row=2, col=1)
    fig.update_xaxes(title_text=xl, row=rows, col=1)
    return fig


def scrub_hist_X_widget(debug, Naxis=1, dt=None, show_slew=True, continuous=True, width=1000):
    """Kernel-backed plotly FigureWidget scrubber -- the one to use for LONG runs.

    Unlike :func:`scrub_hist_X` (which embeds every frame in the cell output and
    bogs down past a few hundred iterations), this keeps the history in Python and
    the slider patches only the *current* frame's data via ``fig.batch_update()``.
    So its cost is O(1) per step regardless of the number of iterations, and it
    stays smooth while you scrub. Returns an ipywidgets VBox -- display it as the
    last expression in a cell.

    Parameters
    ----------
    debug, Naxis, dt, show_slew
        As in :func:`scrub_hist_X`.
    continuous : bool
        If True (default) the plot updates live while you drag. Set False if a
        very large N makes dragging laggy (then it updates on release).
    width : int
        Figure width in px; the slider is stretched to match it.

    Notes
    -----
    Needs ``plotly`` + ``ipywidgets``. FigureWidget rendering in VS Code can be
    finicky; if it shows blank, that's the widget stack, not this code -- say so
    and I'll give you a bqplot/fastplotlib version instead.
    """
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        from ipywidgets import IntSlider, VBox, Layout
    except ImportError as e:
        raise ImportError("scrub_hist_X_widget needs plotly and ipywidgets.") from e

    g, n, N = _reshape_hist_X(debug, Naxis)
    t = np.arange(N) * dt * 1e3 if dt else np.arange(N).astype(float)
    xl = "time [ms]" if dt else "sample"
    slew = show_slew and (dt is not None) and N > 1
    if slew:
        s = np.diff(g, axis=2) / dt
        ts = t[:-1]

    bval = _as1d(debug.get("hist_bvalue"))
    feas = _as1d(debug.get("hist_all_feas"))
    best = int(debug.get("best_feasible_iter", -1))

    gpad = 0.05 * (g.max() - g.min() + 1e-12)
    grange = [g.min() - gpad, g.max() + gpad]
    if slew:
        spad = 0.05 * (s.max() - s.min() + 1e-12)
        srange = [s.min() - spad, s.max() + spad]

    def title(i):
        p = [f"iter {i}/{n - 1}"]
        if bval is not None and i < len(bval):
            p.append(f"b={bval[i]:.1f}")
        if feas is not None and i < len(feas):
            p.append("feasible" if feas[i] else "INFEASIBLE")
        if i == best:
            p.append("&#8592; returned")
        return "   ".join(p)

    rows = 2 if slew else 1
    base = make_subplots(rows=rows, cols=1, shared_xaxes=True, vertical_spacing=0.08)
    for j in range(Naxis):
        base.add_trace(go.Scatter(x=t, y=g[0, j], mode="lines", name=f"axis {j}"), row=1, col=1)
    if slew:
        for j in range(Naxis):
            base.add_trace(go.Scatter(x=ts, y=s[0, j], mode="lines", showlegend=False), row=2, col=1)
    base.update_layout(width=width, height=300 * rows, title_text=title(0), margin=dict(t=60))
    base.update_yaxes(title_text="gradient", range=grange, row=1, col=1)
    if slew:
        base.update_yaxes(title_text="slew", range=srange, row=2, col=1)
    base.update_xaxes(title_text=xl, row=rows, col=1)

    fig = go.FigureWidget(base)
    slider = IntSlider(min=0, max=n - 1, value=0, description="iter", continuous_update=continuous,
                       layout=Layout(width=f"{width}px"), style={"description_width": "30px"})

    def on_change(change):
        i = change["new"]
        with fig.batch_update():  # patch only this frame's data -- nothing is embedded/re-rendered wholesale
            for j in range(Naxis):
                fig.data[j].y = g[i, j]
            if slew:
                for j in range(Naxis):
                    fig.data[Naxis + j].y = s[i, j]
            fig.layout.title.text = title(i)

    slider.observe(on_change, names="value")
    return VBox([slider, fig])


def scrub_hist_X_mpl(debug, Naxis=1, dt=None, show_slew=True):
    """Matplotlib fallback scrubber over hist_X (inline backend + ipywidgets).

    Redraws each frame from scratch on the default inline backend, so it never
    corrupts other figures the way ``%matplotlib widget`` (ipympl) can. Slightly
    less smooth than :func:`scrub_hist_X` (plotly) but has no extra deps beyond
    matplotlib + ipywidgets.
    """
    try:
        import matplotlib.pyplot as plt
        from ipywidgets import interact, IntSlider
    except ImportError as e:
        raise ImportError(
            "scrub_hist_X_mpl needs matplotlib and ipywidgets "
            "(`pixi add --feature dev ipywidgets`)."
        ) from e

    g, n, N = _reshape_hist_X(debug, Naxis)
    t = np.arange(N) * dt * 1e3 if dt else np.arange(N).astype(float)
    xl = "time [ms]" if dt else "sample"
    gmin, gmax = g.min(), g.max()
    gpad = 0.05 * (gmax - gmin + 1e-12)
    slew = show_slew and (dt is not None) and N > 1
    if slew:
        s = np.diff(g, axis=2) / dt
        smin, smax = s.min(), s.max()
        spad = 0.05 * (smax - smin + 1e-12)
        ts = t[:-1]

    bval = _as1d(debug.get("hist_bvalue"))
    feas = _as1d(debug.get("hist_all_feas"))
    best = int(debug.get("best_feasible_iter", -1))

    def plot(i):
        nr = 2 if slew else 1
        fig, ax = plt.subplots(nr, 1, figsize=(9, 3.2 * nr), squeeze=False)
        ax = ax[:, 0]
        for j in range(Naxis):
            ax[0].plot(t, g[i, j], lw=1, label=f"axis {j}")
        ax[0].set_ylim(gmin - gpad, gmax + gpad)
        ax[0].set_ylabel("gradient")
        ax[0].set_xlabel(xl)
        ax[0].grid(alpha=0.3)
        if Naxis > 1:
            ax[0].legend(loc="upper right", fontsize=8)
        if slew:
            for j in range(Naxis):
                ax[1].plot(ts, s[i, j], lw=0.8)
            ax[1].set_ylim(smin - spad, smax + spad)
            ax[1].set_ylabel("slew")
            ax[1].set_xlabel(xl)
            ax[1].grid(alpha=0.3)
        p = [f"iter {i}/{n - 1}"]
        if bval is not None and i < len(bval):
            p.append(f"b={bval[i]:.1f}")
        if feas is not None and i < len(feas):
            p.append("feasible" if feas[i] else "INFEASIBLE")
        if i == best:
            p.append("<-- returned")
        ax[0].set_title("    ".join(p))
        fig.tight_layout()
        plt.show()

    interact(plot, i=IntSlider(min=0, max=n - 1, value=0, description="iter", continuous_update=False))


# --------------------------------------------------------------------------- #
# solve-convergence dashboard
# --------------------------------------------------------------------------- #
def plot_convergence(debug, op_names=None, figsize=(12, 9)):
    """Multi-panel matplotlib dashboard of the solve's convergence history.

    Auto-detects which histories are populated and draws a panel for each:
    b-value vs feasibility, primal/dual residuals, per-operator feasibility
    distance, ADMM weights + gamma_x, inner CG iterations/residual, and the
    objective-vs-constraint pull balance.

    Parameters
    ----------
    debug : dict
        Output of ``solver.get_debug()``.
    op_names : list[str] or None
        Optional operator names for the per-operator panels (order must match
        the solver's operator order). Falls back to ``op0, op1, ...``.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError as e:
        raise ImportError("plot_convergence needs matplotlib.") from e

    bval = _as1d(debug.get("hist_bvalue"))
    feas = _as1d(debug.get("hist_all_feas"))
    best = int(debug.get("best_feasible_iter", -1))
    r_prim = _as2d(debug.get("hist_r_prim"))
    r_dual = _as2d(debug.get("hist_r_dual"))
    r_feas = _as2d(debug.get("hist_r_feas"))
    weight = _as2d(debug.get("hist_weight"))
    gamma_x = _as1d(debug.get("hist_gamma_x"))
    cg_iter = _as1d(debug.get("hist_cg_iter"))
    cg_rnorm = _as1d(debug.get("hist_cg_rnorm"))
    cg_bnorm0 = _as1d(debug.get("hist_cg_bnorm0"))
    obj_pull = _as1d(debug.get("hist_obj_pull"))
    con_pull = _as1d(debug.get("hist_con_pull"))

    # Each panel is (condition, draw_fn). Only those with data are laid out.
    panels = []

    if bval is not None:
        def _p_bval(ax):
            it = np.arange(len(bval))
            ax.plot(it, bval, "-", color="0.6", lw=1, zorder=1)
            if feas is not None:
                m = min(len(it), len(feas))
                fb = feas[:m].astype(bool)
                ax.scatter(it[:m][fb], bval[:m][fb], s=10, c="tab:green", label="feasible", zorder=2)
                ax.scatter(it[:m][~fb], bval[:m][~fb], s=10, c="tab:red", label="infeasible", zorder=2)
                ax.legend(loc="lower right", fontsize=8)
            if 0 <= best < len(bval):
                ax.axvline(best, ls="--", color="k", alpha=0.6)
                ax.annotate("returned", (best, bval[best]), fontsize=8,
                            xytext=(4, 4), textcoords="offset points")
            ax.set_title("b-value vs feasibility")
            ax.set_xlabel("iter"); ax.set_ylabel("b-value")
        panels.append(_p_bval)

    if r_prim is not None or r_dual is not None:
        def _p_resid(ax):
            if r_prim is not None:
                ax.semilogy(np.nanmax(r_prim, axis=1), label="max primal", color="tab:blue")
            if r_dual is not None:
                ax.semilogy(np.nanmax(r_dual, axis=1), label="max dual", color="tab:orange")
            ax.set_title("ADMM residuals (max over ops)")
            ax.set_xlabel("iter"); ax.set_ylabel("residual"); ax.legend(fontsize=8)
        panels.append(_p_resid)

    if r_feas is not None:
        def _p_feas(ax):
            for i in range(r_feas.shape[1]):
                ax.semilogy(np.clip(r_feas[:, i], 1e-16, None), lw=1, label=_op_label(i, op_names))
            ax.set_title("per-op feasibility distance")
            ax.set_xlabel("iter"); ax.set_ylabel("r_feas")
            if r_feas.shape[1] <= 10:
                ax.legend(fontsize=7, ncol=2)
        panels.append(_p_feas)

    if weight is not None:
        def _p_weight(ax):
            for i in range(weight.shape[1]):
                ax.semilogy(weight[:, i], lw=1, label=_op_label(i, op_names))
            ax.set_title("operator weights (rho)")
            ax.set_xlabel("iter"); ax.set_ylabel("weight")
            if weight.shape[1] <= 10:
                ax.legend(fontsize=7, ncol=2)
            if gamma_x is not None:
                ax2 = ax.twinx()
                ax2.plot(gamma_x, color="0.5", ls=":", lw=1)
                ax2.set_ylabel("gamma_x", color="0.5")
        panels.append(_p_weight)

    if cg_iter is not None or cg_rnorm is not None:
        def _p_cg(ax):
            if cg_iter is not None:
                ax.plot(cg_iter, color="tab:purple", lw=1, label="CG iters")
                ax.set_ylabel("CG iters", color="tab:purple")
            ax.set_title("inner CG solve"); ax.set_xlabel("iter")
            if cg_rnorm is not None:
                ax2 = ax.twinx()
                rn = cg_rnorm.copy()
                if cg_bnorm0 is not None:
                    m = min(len(rn), len(cg_bnorm0))
                    denom = np.where(cg_bnorm0[:m] > 0, cg_bnorm0[:m], 1.0)
                    rn = rn[:m] / denom
                    lbl = "CG rnorm / bnorm0"
                else:
                    lbl = "CG rnorm"
                ax2.semilogy(np.clip(rn, 1e-20, None), color="tab:gray", lw=1)
                ax2.set_ylabel(lbl, color="tab:gray")
        panels.append(_p_cg)

    if obj_pull is not None or con_pull is not None:
        def _p_pull(ax):
            if obj_pull is not None:
                ax.semilogy(np.clip(obj_pull, 1e-20, None), label="obj pull", color="tab:green")
            if con_pull is not None:
                ax.semilogy(np.clip(con_pull, 1e-20, None), label="con pull", color="tab:red")
            ax.set_title("objective vs constraint pull")
            ax.set_xlabel("iter"); ax.set_ylabel("||pull||"); ax.legend(fontsize=8)
        panels.append(_p_pull)

    if not panels:
        print("No plottable histories found -- set solver.extra_debug = True before solve().")
        return None

    ncols = 2
    nrows = (len(panels) + ncols - 1) // ncols
    fig, axarr = plt.subplots(nrows, ncols, squeeze=False, figsize=figsize, layout="tight")
    axes = axarr.flatten()
    for ax, draw in zip(axes, panels):
        draw(ax)
        ax.grid(alpha=0.3)
        ax.spines["top"].set_visible(False)
    for ax in axes[len(panels):]:
        ax.set_visible(False)
    plt.show()
    return fig
