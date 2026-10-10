"""Home-page figures, in a light and a dark version (the page shows the one matching the theme).

    pixi run -e dev python docs/img/make_figures.py

hero_*.png: the same b-value at a shorter TE (conventional symmetric timing vs GrOpt).
quickstart_*.png: the waveform from the home-page quick start.
"""
from dataclasses import replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import gropt
from gropt.diffusion import DiffParams, te_search, to_raster

OUT = Path(__file__).parent
CFG = DiffParams(TE=80e-3, T_90=3e-3, T_180=5e-3, T_readout=16e-3, dt=100e-6, gmax=0.08, smax=100.0, MMT=1)
TARGET_B = 1000.0

# highlight GrOpt, gray the comparison (validated: CVD and contrast pass in both modes)
THEMES = {
    "light": dict(grop="#9b1b1b", conv="#80858f", text="#2b2b2b", muted="#6b6f76", block="#e9eaee"),
    "dark": dict(grop="#e0524c", conv="#858b94", text="#e8e8e8", muted="#a3a8b0", block="#2f333d"),
}


def style(ax, c):
    ax.patch.set_alpha(0.0)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(c["muted"])
    ax.tick_params(colors=c["muted"], labelsize=9.5)


def search(mode):
    out = te_search(replace(CFG, diff_mode=mode), target_b=TARGET_B)
    hw = to_raster(out["result"], 10e-6)
    return out["TE"], out["result"]["bvalue"], np.asarray(hw["X"])


def draw_hero(rows, c, path):
    te_max = max(te for _, te, _, _, _ in rows)
    fig, axs = plt.subplots(2, 1, figsize=(10, 3.9), sharex=True, gridspec_kw={"hspace": 0.35})
    fig.patch.set_alpha(0.0)
    for ax, (label, te, b, g, color) in zip(axs, rows):
        style(ax, c)
        t = np.arange(g.size) * 10e-3                       # ms from the excitation
        gm = g * 1e3                                         # mT/m
        # RF pulses and readout, as quiet blocks behind the waveform
        for t0, t1, name in ((0.0, CFG.T_90 * 1e3, "90°"),
                             ((te - CFG.T_180) * 500, (te + CFG.T_180) * 500, "180°"),
                             ((te - CFG.T_readout) * 1e3, te * 1e3 + 4, "readout")):
            ax.axvspan(t0, t1, color=c["block"], zorder=0, lw=0)
            ax.text((t0 + t1) / 2, 92, name, ha="center", va="bottom", fontsize=9, color=c["muted"])
        ax.axhline(0, color=c["muted"], lw=0.6, zorder=1)
        ax.fill_between(t, gm, color=color, alpha=0.18, lw=0, zorder=2)
        ax.plot(t, gm, color=color, lw=2, zorder=3, solid_capstyle="round")
        ax.axvline(te * 1e3, color=color, lw=1.5, ls=(0, (4, 3)), zorder=3)
        ax.text(te * 1e3 + 1.0, -78, f"TE {te * 1e3:.1f} ms", color=c["text"], fontsize=10.5, va="bottom",
                fontweight="bold")
        title = ax.text(0, 118, label, color=c["text"], fontsize=12, fontweight="bold", va="bottom")
        ax.annotate(f"   b = {b:.0f} s/mm²", xy=(1, 0), xycoords=title, va="bottom", color=c["muted"],
                    fontsize=10.5)
        ax.set_ylim(-95, 112)
        ax.set_yticks([])
    # on the GrOpt row: where the conventional TE was, and the difference
    (_, te_c, *_), (_, te_g, *_) = rows
    ax = axs[-1]
    ax.axvline(te_c * 1e3, color=c["conv"], lw=1.2, ls=(0, (4, 3)), alpha=0.6, zorder=3)
    ax.annotate("", xy=(te_g * 1e3, 55), xytext=(te_c * 1e3, 55),
                arrowprops=dict(arrowstyle="-|>", color=c["text"], lw=1.3, shrinkA=0, shrinkB=0))
    ax.text((te_g + te_c) * 500, 62, f"{(te_c - te_g) * 1e3:.1f} ms shorter", ha="center", va="bottom",
            color=c["text"], fontsize=10.5, fontweight="bold")
    axs[-1].set_xlabel("time from excitation [ms]", color=c["muted"], fontsize=10)
    axs[-1].set_xlim(-1, te_max * 1e3 + 6)
    fig.savefig(path, dpi=200, bbox_inches="tight", transparent=True)
    plt.close(fig)


def quickstart():
    # the home-page quick start, unchanged
    gparams = gropt.GroptParams()
    gparams.dt = 10e-6
    gparams.N = 112
    gparams.vec_init_simple()
    gparams.add_gmax(.08)
    gparams.add_smax(200)
    gparams.add_moment(0, 14.0)
    gparams.add_moment(1, 0.0)
    return gropt.solve(gparams).X, gparams.dt


def draw_quickstart(g, dt, c, path):
    fig, ax = plt.subplots(figsize=(5.2, 3.0))
    fig.patch.set_alpha(0.0)
    style(ax, c)
    t = np.arange(g.size) * dt * 1e3
    for lim in (-80, 80):
        ax.axhline(lim, color=c["muted"], lw=0.8, ls=(0, (2, 3)))
    ax.text(t[-1], 80, "gmax", ha="right", va="bottom", fontsize=9, color=c["muted"])
    ax.axhline(0, color=c["muted"], lw=0.6)
    ax.fill_between(t, g * 1e3, color=c["grop"], alpha=0.18, lw=0)
    ax.plot(t, g * 1e3, color=c["grop"], lw=2, solid_capstyle="round")
    ax.set_xlabel("t [ms]", color=c["muted"], fontsize=10)
    ax.set_ylabel("G [mT/m]", color=c["muted"], fontsize=10)
    ax.spines["left"].set_visible(True)
    ax.spines["left"].set_color(c["muted"])
    ax.text(0.98, 0.72, "M0 = 14 mT·ms/m\nM1 = 0", transform=ax.transAxes, ha="right", va="center",
            fontsize=9.5, color=c["text"])
    fig.savefig(path, dpi=200, bbox_inches="tight", transparent=True)
    plt.close(fig)


if __name__ == "__main__":
    plt.rcParams.update({"font.family": ["Arial", "DejaVu Sans"], "font.size": 11})
    te_c, b_c, g_c = search("conventional")
    te_g, b_g, g_g = search("gropt")
    print(f"conventional TE {te_c * 1e3:.2f} ms b {b_c:.1f}; GrOpt TE {te_g * 1e3:.2f} ms b {b_g:.1f}")
    g_q, dt_q = quickstart()
    for name, c in THEMES.items():
        rows = [("Conventional (symmetric)", te_c, b_c, g_c, c["conv"]), ("GrOpt", te_g, b_g, g_g, c["grop"])]
        draw_hero(rows, c, OUT / f"hero_{name}.png")
        draw_quickstart(g_q, dt_q, c, OUT / f"quickstart_{name}.png")
