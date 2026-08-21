import json
import pickle
from dataclasses import fields, replace
from pathlib import Path

import numpy as np

# ── DiffParams fields that DEFINE the problem (the optimum) -> NOT part of a solve recipe. Everything
#    else in DiffParams is a solve knob and is captured automatically, so new knobs need no edits here.
PROBLEM_FIELDS = {
    "TE", "T_90", "T_180", "T_readout", "T_pre", "dt", "diff_mode",                      # layout
    "gmax", "smax",                                                                      # hardware limits
    "MMT", "moment_tol", "pns_lim", "cns_lim", "safe_params", "concomitant", "eddy_lam", # constraints
    "bvalue", "bval_mode", "bval_min",                                                   # objective / target
    "jerk_lam", "basin_same_sign", "basin_window", "basin_eps",                          # shape the waveform
}

def _recipe_fields(cfg):
    return [f.name for f in fields(cfg) if f.name not in PROBLEM_FIELDS]

def _native(o):   # random_cfgs emits np.float64/np.int64 -> plain python for JSON + clean round-trip
    if isinstance(o, np.integer):  return int(o)
    if isinstance(o, np.floating): return float(o)
    if isinstance(o, np.bool_):    return bool(o)
    raise TypeError(f"not JSON-serializable: {type(o)}")


# ── RECIPE = static "how to solve" config. JSON, portable/cross-language, reusable across problems. ───
def save_recipe(path, name, cfg, scfg, *, description=""):
    """Add or overwrite a named recipe in a JSON library file.

    Stores solve knobs only: all of ``SolverCfg`` plus the non-problem ``DiffParams``
    fields. The geometry and the warmstart are not saved.

    Parameters
    ----------
    path : str or pathlib.Path
        JSON library file to create or update.
    name : str
        Key under which the recipe is stored; overwrites any existing entry.
    cfg : DiffParams
        Problem config; only its non-problem (solve-knob) fields are saved.
    scfg : SolverCfg
        Solver config; every field is saved.
    description : str, optional
        Human-readable note stored alongside the recipe.

    Returns
    -------
    pathlib.Path
        The library file that was written.
    """
    entry = {
        "description": description,
        "diff":   {f: getattr(cfg, f) for f in _recipe_fields(cfg)},
        "solver": {f.name: getattr(scfg, f.name) for f in fields(scfg)},
    }
    p = Path(path)
    lib = json.loads(p.read_text()) if p.exists() else {}
    lib[name] = entry
    p.write_text(json.dumps(lib, indent=2, sort_keys=True, default=_native))
    return p

def load_recipe(path, name, base_cfg, base_scfg):
    """Overlay a recipe's tuning onto your problem's base configs.

    The problem geometry and hardware limits (``TE``, ``dt``, ``gmax``, ...) in the base
    configs are kept; only the recipe's solve knobs are applied on top.

    Parameters
    ----------
    path : str or pathlib.Path
        JSON library file to read.
    name : str
        Recipe key to load.
    base_cfg : DiffParams
        Base problem config to overlay the recipe's ``diff`` fields onto.
    base_scfg : SolverCfg
        Base solver config to overlay the recipe's ``solver`` fields onto.

    Returns
    -------
    cfg : DiffParams
        ``base_cfg`` with the recipe's solve knobs applied.
    scfg : SolverCfg
        ``base_scfg`` with the recipe's solver settings applied.
    """
    r = json.loads(Path(path).read_text())[name]
    return replace(base_cfg, **r["diff"]), replace(base_scfg, **r["solver"])

def list_recipes(path):
    """Return a mapping of recipe name to description for a library file.

    Parameters
    ----------
    path : str or pathlib.Path
        JSON library file to read.

    Returns
    -------
    dict
        ``{name: description}`` for every recipe in the file.
    """
    return {k: v.get("description", "") for k, v in json.loads(Path(path).read_text()).items()}


# ── WARMSTART = dynamic solve state (a converged snapshot to hot-start a nearby problem). Separate
#    binary file, per-solve, pairs with a recipe of the same name. ──────────────────────────────────
def save_warmstart(folder, name, result):
    """Save one solve's converged warmstart snapshot as its own file.

    The snapshot (duals, primal, and adapted weights) is pickled to
    ``{folder}/{name}.warmstart.pkl``.

    Parameters
    ----------
    folder : str or pathlib.Path
        Directory to write into; created if it does not exist.
    name : str
        Base name for the snapshot file; pairs with a recipe of the same name.
    result : dict
        A full solve result; must contain a ``"warmstart"`` snapshot.

    Returns
    -------
    pathlib.Path
        The snapshot file that was written.

    Raises
    ------
    ValueError
        If ``result`` has no ``"warmstart"`` (e.g. a minimized result).
    """
    ws = result.get("warmstart")
    if ws is None:
        raise ValueError("result has no 'warmstart' (minimized result?) -- keep the full result to save it")
    d = Path(folder); d.mkdir(parents=True, exist_ok=True)
    p = d / f"{name}.warmstart.pkl"
    p.write_bytes(pickle.dumps(ws))
    return p

def load_warmstart(folder, name):
    """Load a saved snapshot to pass as ``solve(..., warmstart=...)``.

    Parameters
    ----------
    folder : str or pathlib.Path
        Directory the snapshot was saved in.
    name : str
        Base name used when the snapshot was saved.

    Returns
    -------
    dict
        The warmstart snapshot.
    """
    return pickle.loads((Path(folder) / f"{name}.warmstart.pkl").read_bytes())