"""Tests for the shipped diffusion recipe library and the recipe API.

    pixi run -e dev test-py
"""
import json
import tempfile
import warnings
from dataclasses import fields, replace
from pathlib import Path

import numpy as np

import gropt
import gropt.diffusion as gd
from gropt.diffusion_recipes import LIBRARY, PROBLEM_FIELDS, RECIPE_FIELDS, load_recipe, save_recipe

LIB = json.loads(LIBRARY.read_text())
PNS = gd.DiffParams(TE=80e-3, T_90=3e-3, T_180=5e-3, T_readout=16e-3, dt=200e-6, gmax=0.08, MMT=1,
                    pns_lim=0.8, cns_lim=0.8)
CONC = gd.DiffParams(TE=60e-3, T_90=3e-3, T_180=5e-3, T_readout=16e-3, dt=400e-6, gmax=0.19, MMT=2,
                     concomitant=True, concomitant_tol=0.01)
# b per role on PNS / CONC, recorded 2026-09-27; chaotic in the last digits, so compared at 2%
REF_B = {"default": (1606.98, 525.09), "best_single": (1609.11, 525.12), "fast": (1608.21, 525.22),
         "fastest": (1603.87, 524.83), "soft_concomitant": (1609.22, 527.34)}


def test_every_field_classified():
    names = {f.name for f in fields(gd.DiffParams)}
    assert not PROBLEM_FIELDS & RECIPE_FIELDS, PROBLEM_FIELDS & RECIPE_FIELDS
    assert names == PROBLEM_FIELDS | RECIPE_FIELDS, (
        f"unclassified: {sorted(names - PROBLEM_FIELDS - RECIPE_FIELDS)}; "
        f"stale: {sorted((PROBLEM_FIELDS | RECIPE_FIELDS) - names)}")


def test_library_complete():
    solver_names = {f.name for f in fields(gd.SolverCfg)}
    for name, e in LIB["recipes"].items():
        assert set(e["diff"]) == RECIPE_FIELDS, (name, set(e["diff"]) ^ RECIPE_FIELDS)
        assert set(e["solver"]) == solver_names, (name, set(e["solver"]) ^ solver_names)
    for role, key in LIB["roles"].items():
        assert key in LIB["recipes"], (role, key)
    for p, v in LIB["portfolios"].items():
        assert all(LIB["roles"].get(m, m) in LIB["recipes"] for m in v["members"]), p
    assert "default" in LIB["roles"]


def test_defaults_are_default_recipe():
    e = LIB["recipes"][LIB["roles"]["default"]]
    dflt = {f.name: f.default for f in fields(gd.DiffParams)}
    assert {k: dflt[k] for k in RECIPE_FIELDS} == e["diff"], "DiffParams defaults != the 'default' recipe"
    assert {f.name: f.default for f in fields(gd.SolverCfg)} == e["solver"], "SolverCfg defaults != 'default'"
    cfg, scfg = gd.apply_recipe(PNS, "default")
    assert cfg == PNS and scfg == gd.SolverCfg()


def test_apply_keeps_problem_and_round_trips():
    cfg, scfg = gd.apply_recipe(replace(CONC, concomitant_project=False), "fast")
    assert cfg.concomitant_project is False and cfg.TE == CONC.TE     # problem fields untouched
    with tempfile.TemporaryDirectory() as d:
        path = Path(d) / "mine.json"
        save_recipe(path, "mine", cfg, scfg)
        entry = json.loads(path.read_text())["mine"]
        assert set(entry["diff"]) == RECIPE_FIELDS
        assert load_recipe(path, "mine", CONC, gd.SolverCfg()) == (replace(cfg, concomitant_project=True), scfg)
        entry["diff"]["concomitant_project"] = False    # older files stored this problem field: skip it
        path.write_text(json.dumps({"old": entry}))
        assert load_recipe(path, "old", CONC, gd.SolverCfg())[0].concomitant_project is True
        del entry["diff"]["safe_lifted"]                # saved before the lifting existed: runs unlifted
        path.write_text(json.dumps({"pre_lift": entry}))
        assert load_recipe(path, "pre_lift", CONC, gd.SolverCfg())[0].safe_lifted is False
        entry["diff"]["no_such_knob"] = 1.0
        path.write_text(json.dumps({"bad": entry}))
        try:
            load_recipe(path, "bad", CONC, gd.SolverCfg())
            raise AssertionError("unknown field accepted")
        except ValueError:
            pass


def test_errors_and_warning():
    for bad, exc in ((lambda: gd.solve(PNS, gd.SolverCfg(), recipe="fast"), ValueError),
                     (lambda: gd.apply_recipe(PNS, "best"), ValueError),
                     (lambda: gd.apply_recipe(PNS, "no_such_recipe"), KeyError)):
        try:
            bad()
            raise AssertionError("no error")
        except exc:
            pass
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        gd.apply_recipe(replace(PNS, w_pns=5.0), "fast")
    assert any("w_pns" in str(x.message) for x in w)


def test_iter_recipes():
    names = [n for n, _, _ in gd.iter_recipes(PNS, "best")]
    assert names == LIB["portfolios"]["best"]["members"], names
    assert [n for n, _, _ in gd.iter_recipes(PNS)] == list(LIB["roles"])
    (name, cfg, scfg), = gd.iter_recipes(PNS, "fast")
    assert (cfg, scfg) == gd.apply_recipe(PNS, "fast")


def test_verify():
    r = gd.solve(PNS)
    v = gd.verify(r, PNS)
    assert r["converged"] and v["feasible"] and 0.9 < v["max_ratio"] <= 1.02, v
    bad = gd.verify({**r, "X": 1.5 * np.asarray(r["X"])}, PNS)
    assert not bad["feasible"] and bad["g_ratio"] > 1.4, bad
    assert not gd.verify({"X": None}, PNS)["feasible"]


def test_portfolio():
    r = gd.solve(CONC, recipe="best", parallel=False)
    rows = r["portfolio"]
    assert len(rows) == len(LIB["portfolios"]["best"]["members"])
    ok = [x["bvalue"] for x in rows if x["converged"] and x["feasible"]]
    assert ok and r["bvalue"] == max(ok) and r["verify"]["feasible"]


def test_regression():
    for role, refs in REF_B.items():
        for cfg, b0 in zip((PNS, CONC), refs, strict=True):
            if role == "soft_concomitant" and cfg.concomitant:
                cfg = replace(cfg, concomitant_project=False)
            r = gd.solve(cfg, recipe=role)
            assert r["converged"] and gd.verify(r, cfg)["feasible"], role
            assert abs(r["bvalue"] - b0) <= 0.02 * b0, (role, r["bvalue"], b0)


if __name__ == "__main__":
    gropt.set_log_level(6)
    tests = [(k, v) for k, v in dict(globals()).items() if k.startswith("test_")]
    for name, fn in tests:
        fn()
        print(f"ok  {name}")
    print(f"{len(tests)} passed")
