"""Diffusion solve recipes and warm-start snapshots.

A recipe is the solve configuration: every ``SolverCfg`` field plus the ``DiffParams`` fields in
``RECIPE_FIELDS``; the problem (``PROBLEM_FIELDS``) always comes from the caller. gropt ships a tuned library
(``diffusion_recipes.json`` next to this file, see :func:`available_recipes`); :func:`save_recipe` writes your
own. A warmstart is the solver snapshot from one solve, pickled to its own ``{name}.warmstart.pkl`` file.
"""

import json
import pickle
import warnings
from dataclasses import fields, replace
from pathlib import Path

import numpy as np

# DiffParams fields that define the problem; never stored in or taken from a recipe.
PROBLEM_FIELDS = {
    "TE", "T_90", "T_180", "T_readout", "T_pre", "dt", "diff_mode",                      # layout
    "gmax", "smax",                                                                      # hardware limits
    "MMT", "moment_tol", "pns_lim", "cns_lim", "safe_params", "concomitant", "eddy_lam", # constraints
    "concomitant_tol", "concomitant_project", "safe_test_axes",                          # constraints
    "bvalue", "bval_mode", "bval_min",                                                   # objective / target
    "jerk_lam", "basin_same_sign", "basin_window", "basin_eps",                          # shape the waveform
    "safe_alpha_exact", "concomitant_exact_quad", "moment_pwl_quad", "bval_pwl_quad",    # model accuracy
}
# DiffParams fields a recipe carries (solve knobs). Every DiffParams field is in exactly one of the two sets.
RECIPE_FIELDS = {
    "w_gmax", "w_smax", "w_moment", "w_bval", "w_pns", "w_cns", "w_slack", "w_concomitant",   # weights
    "w_eddy", "w_jerk", "moment_project", "eddy_project", "bval_obj_weight", "bval_max_scale",
    "safe_eps", "safe_signed13", "safe_lifted",                                          # SAFE formulation
    "x0_mode", "x0_amp", "x0_invert", "x0_periods", "x0_project",                        # initial waveform
}
LIBRARY = Path(__file__).with_name("diffusion_recipes.json")
CURATED = 2   # "schema" of a curated library (roles, portfolios, dated entries)

def _recipe_fields(cfg):
    return [f.name for f in fields(cfg) if f.name in RECIPE_FIELDS]

def _native(o):   # numpy scalars -> plain Python for JSON
    if isinstance(o, np.integer):  return int(o)
    if isinstance(o, np.floating): return float(o)
    if isinstance(o, np.bool_):    return bool(o)
    raise TypeError(f"not JSON-serializable: {type(o)}")

def _read_library(path=None):
    """Read a library file as ``{"roles", "portfolios", "recipes"}``; a save_recipe file has no roles."""
    lib = json.loads(Path(path or LIBRARY).read_text())
    if lib.get("schema") == CURATED:
        return lib
    return {"roles": {}, "portfolios": {}, "recipes": lib}

def _entry(lib, name, path=None):
    """Return the entry a role or entry name points at; portfolio names raise (they run through solve)."""
    if name in lib["portfolios"]:
        msg = f"{name!r} is a portfolio; use gropt.diffusion.solve(cfg, recipe={name!r})"
        raise ValueError(msg)
    key = lib["roles"].get(name, name)
    if key not in lib["recipes"]:
        have = sorted({*lib["roles"], *lib["portfolios"]} or lib["recipes"])
        msg = f"no recipe {name!r} in {path or LIBRARY.name} (have: {', '.join(have)})"
        raise KeyError(msg)
    return lib["recipes"][key]

def _overlay(cfg, scfg, entry, name, warn):
    """``cfg`` / ``scfg`` with a library entry's knobs applied (problem fields in older files are skipped)."""
    diff = {"safe_lifted": False,   # a recipe saved before the lifting existed was tuned unlifted
            **{k: v for k, v in entry.get("diff", {}).items() if k not in PROBLEM_FIELDS}}
    solver = entry.get("solver", {})
    unknown = sorted((set(diff) - RECIPE_FIELDS) | (set(solver) - {f.name for f in fields(scfg)}))
    if unknown:
        msg = f"recipe {name!r} has unknown fields {unknown} (saved by a newer gropt?)"
        raise ValueError(msg)
    if warn:
        default = {f.name: f.default for f in fields(cfg)}
        for k, v in diff.items():
            cur = getattr(cfg, k)
            if cur != default[k] and cur != v:
                warnings.warn(f"recipe {name!r} overrides {k}={cur!r} with {v!r}", stacklevel=3)
    return replace(cfg, **diff), replace(scfg, **solver)


# Recipes: solve settings in a JSON library, reusable across problems.
def available_recipes(all=False, library=None):  # noqa: A002
    """List the recipes you can pass to ``solve(cfg, recipe=name)``: what each is for and how it scored.

    Scores come from the tuning grid (PNS/CNS on and off, concomitant off / exact / band, 100-400 us, plus a
    held-out grid): ``score`` is the mean over problems of b / best known b, ``worst`` the weakest use case,
    ``time`` the total solve time relative to ``default``. A portfolio runs several recipes and keeps the
    best verified result (see :func:`gropt.diffusion.solve`).

    Parameters
    ----------
    all : bool, optional
        Also list every dated entry, including portfolio members that are not a named role.
    library : str or pathlib.Path, optional
        A library file; default the one shipped with gropt.

    Returns
    -------
    RecipeMenu
        ``{name: info dict}``; prints as a table.
    """
    def row(kind, entry, e):
        return {"kind": kind, "entry": entry, "description": e.get("description", ""), **e.get("scores", {})}

    lib = _read_library(library)
    menu = RecipeMenu({r: row("recipe", k, lib["recipes"][k]) for r, k in lib["roles"].items()})
    menu.update({n: row(f"portfolio of {len(p['members'])}", ", ".join(p["members"]), p)
                 for n, p in lib["portfolios"].items()})
    if all or not menu:
        for k, e in lib["recipes"].items():
            menu.setdefault(k, row("recipe", k, e))
    return menu

class RecipeMenu(dict):
    """``available_recipes()`` result: a dict that prints as a table."""

    def __repr__(self):
        """Plain-text table."""
        cols = (("kind", "{}"), ("score", "{:.3f}"), ("worst", "{:.3f}"), ("time", "{:.2f}x"))
        rows = [["name", *(c for c, _ in cols), "what for"]]
        rows += [[n, *(f.format(r[c]) if c in r else "-" for c, f in cols), r["description"]]
                 for n, r in self.items()]
        w = [max(len(r[i]) for r in rows) for i in range(len(cols) + 1)]
        return "\n".join("  ".join(c.ljust(n) for c, n in zip(r, w, strict=False)) + "  " + r[-1]
                         for r in rows)

def apply_recipe(cfg, name, *, library=None, warn=True):
    """Your problem with a named recipe's solve knobs applied.

    Parameters
    ----------
    cfg : DiffParams
        The problem; its problem fields (timing, limits, constraints) are kept.
    name : str
        A role (``"default"``, ``"fast"``, ...) or dated entry from :func:`available_recipes`.
    library : str or pathlib.Path, optional
        A library file; default the one shipped with gropt.
    warn : bool, optional
        Warn when the recipe replaces a knob you changed from its default. Default True.

    Returns
    -------
    cfg : DiffParams
        ``cfg`` with the recipe's knobs.
    scfg : SolverCfg
        The recipe's solver settings.

    Raises
    ------
    KeyError
        Unknown name.
    ValueError
        ``name`` is a portfolio (pass it to ``solve(cfg, recipe=name)``), or the entry has fields this gropt
        does not know.
    """
    from gropt.diffusion import SolverCfg  # noqa: PLC0415 -- gropt.diffusion imports this module

    return _overlay(cfg, SolverCfg(), _entry(_read_library(library), name, library), name, warn)

def iter_recipes(cfg, names=None, *, library=None):
    """Yield your problem under each of several recipes, for your own loops (e.g. a TE search per recipe).

    ``solve(cfg, recipe="best")`` already runs a portfolio and keeps the best verified result; iterate when
    you need something else per recipe::

        for name, c, s in iter_recipes(cfg, "best"):
            out = te_search(c, s, target_b=1000.0)

    Parameters
    ----------
    cfg : DiffParams
        The problem.
    names : str or list of str, optional
        A portfolio (its members, in order), one recipe, or a list of recipes; default every named recipe
        (the roles of :func:`available_recipes`).
    library : str or pathlib.Path, optional
        A library file; default the one shipped with gropt.

    Yields
    ------
    name : str
        The recipe (a portfolio yields its member entries).
    cfg : DiffParams
        ``cfg`` with the recipe's knobs.
    scfg : SolverCfg
        The recipe's solver settings.
    """
    lib = _read_library(library)
    if names is None:
        names = list(lib["roles"]) or list(lib["recipes"])
    elif isinstance(names, str):
        names = lib["portfolios"][names]["members"] if names in lib["portfolios"] else [names]
    for i, name in enumerate(names):
        yield (name, *apply_recipe(cfg, name, library=library, warn=(i == 0)))

def save_recipe(path, name, cfg, scfg, *, description=""):
    """Add or overwrite a named recipe in a JSON library file.

    Stores solve knobs only: all of ``SolverCfg`` plus the ``RECIPE_FIELDS`` of ``DiffParams``. The
    problem and the warmstart are not saved.

    Parameters
    ----------
    path : str or pathlib.Path
        JSON library file to create or update.
    name : str
        Key under which the recipe is stored; overwrites any existing entry.
    cfg : DiffParams
        Problem config; only its solve-knob fields are saved.
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
    if lib.get("schema") == CURATED:
        msg = f"{p} is a curated library; save to your own file"
        raise ValueError(msg)
    lib[name] = entry
    p.write_text(json.dumps(lib, indent=2, sort_keys=True, default=_native))
    return p

def load_recipe(path, name, base_cfg, base_scfg):
    """Overlay a recipe's tuning onto your problem's base configs.

    The problem fields in the base configs (``TE``, ``dt``, ``gmax``, ...) are kept; only the recipe's
    solve knobs are applied on top. For the shipped library, :func:`apply_recipe` is shorter.

    Parameters
    ----------
    path : str or pathlib.Path
        JSON library file to read (your own, or ``LIBRARY``).
    name : str
        Recipe key (or, in the shipped library, a role).
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
    return _overlay(base_cfg, base_scfg, _entry(_read_library(path), name, path), name, warn=False)

def list_recipes(path=None):
    """Return a mapping of recipe name to description for a library file.

    Parameters
    ----------
    path : str or pathlib.Path, optional
        JSON library file to read; default the shipped library (roles and portfolios; see
        :func:`available_recipes` for scores).

    Returns
    -------
    dict
        ``{name: description}``.
    """
    return {k: v["description"] for k, v in available_recipes(library=path).items()}


# Warmstarts: per-solve snapshots (pickle) for hot-starting a nearby problem.
def save_warmstart(folder, name, result):
    """Pickle one solve's warmstart snapshot to ``{folder}/{name}.warmstart.pkl``.

    The snapshot holds the primal, duals, and adapted weights.

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
        If ``result`` has no ``"warmstart"`` (e.g. a failed solve or a stripped result).
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