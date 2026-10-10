# Diffusion

`gropt.diffusion` builds and solves diffusion encoding waveforms from a `DiffParams` problem description.

```python
import gropt.diffusion as gd

cfg = gd.DiffParams(TE=80e-3, T_90=3e-3, T_180=5e-3, T_readout=16e-3, gmax=0.08, smax=200.0,
                    MMT=1, pns_lim=0.8, cns_lim=0.8)
res = gd.solve(cfg)              # the tuned default recipe
res["bvalue"], res["X"]
```

## Recipes

A recipe is a set of solve settings (constraint weights, initial waveform, solver knobs) that has been tuned
over many problems: PNS/CNS on and off, concomitant compensation off, exact or within a band, and a range of
timings and hardware limits. The problem itself (timing, limits, constraints) always comes from your
`DiffParams`. The `DiffParams` and `SolverCfg` defaults already are the `default` recipe, so a recipe only
needs naming to trade speed against b-value:

```python
gd.available_recipes()                      # the menu: what each is for, score, worst case, time
res = gd.solve(cfg, recipe="fast")          # about 98.5% of the best known b in ~0.3x the time
res = gd.solve(cfg, recipe="best")          # a portfolio: 3 recipes, keeps the best verified result
cfg2, scfg2 = gd.apply_recipe(cfg, "fast")  # the configs, to adjust with dataclasses.replace
for name, c, s in gd.iter_recipes(cfg, "best"):   # your own loop over a portfolio's recipes
    out = gd.te_search(c, s, target_b=1000.0)
```

`score` in the menu is the mean over the tuning problems of b / best known b; `worst` is the weakest use case;
`time` is relative to `default`. Named recipes (roles) point at dated entries that never change, so after a
re-tune `recipe="univ_20260927_a"` still reproduces an old result (`available_recipes(all=True)` lists them).

**Portfolios** run several complementary recipes, in parallel by default, and return the highest b that passes
`verify`. They escape the b-value basins a single recipe can fall into. The result adds `recipe` (the winner),
`verify` and `portfolio` (b, feasibility and time for every member). Use them for a final waveform, not inside
a TE search.

**Concomitant compensation**: `concomitant_project=True` (the default) enforces an exact balance;
`concomitant_project=False` with `concomitant_tol` allows a band. Both are problem settings, and every recipe
except `soft_concomitant` is tuned to work in both modes.

**Your own recipes**: `save_recipe(path, name, cfg, scfg)` stores the solve knobs of a configuration you tuned,
and `solve(cfg, recipe=name, library=path)` or `load_recipe` uses it. The C++ standalone (`gropt/src`) reads the
same files, including the shipped `gropt/diffusion_recipes.json` (portfolios are Python only).

## Checking a waveform

`converged` means the solver's own constraint checks passed on the solve raster. `verify(result, cfg)` re-checks
every limit independently from the waveform, as the scanner plays it (PNS/CNS on the 10 us hardware raster, exact
concomitant integral), and returns the worst ratio per limit and `feasible`. It guards against a coarse-raster
model that under-reports (`safe_alpha_exact=False`) or a faulty check; `te_search` and portfolios require it.

## Minimum TE and sweeps

`te_search(cfg, target_b=...)` finds the shortest TE whose waveform reaches the target b-value by bisection over
TE, using a single recipe (`recipe="fast"` roughly halves the search time). TE is quantized to `dt`, so the search
runs over the `dt` bins and solves each in its middle; the returned `TE` and `cfg` are the ones `result` was solved
for, so `cfg` re-solves to the same waveform. It uses the exact piecewise-linear b-value (`bval_pwl_quad=True`), so
`b >= target_b` also holds once `to_raster` moves the waveform to the hardware raster.

`sweep(cfg, cfg_grid={...}, scfg_grid={...})` solves every combination of the listed `DiffParams` and `SolverCfg`
values, in parallel processes by default; `sweep_points` takes an explicit list of points instead.

## Solve raster and hardware raster

Solve on a coarse raster (100 to 400 us) and move the result to the 10 us hardware raster with `to_raster`. The
resampling is exact, since the scanner plays the linear interpolant; the cost of a coarse raster is that the RF
gaps are rounded outward to its grid (about 1% of b at 200 us, 5% at 400 us). The solver's own `bvalue` uses a
rectangle rule that reads slightly high on a coarse raster; `bval_pwl_quad=True` makes it exact.

## Under the hood

`build_gparams(cfg)` returns the `GroptParams` that `solve` builds, and `make_solver(scfg)` the configured
`SolverGroptSDMM`, so you can add your own constraints and keep the tuned settings. See the
[examples](../examples/diffusion/diffusion.ipynb), including the
[low-level API](../examples/diffusion/diffusion_constraint.ipynb).

## API

::: gropt.diffusion.solve

::: gropt.diffusion.available_recipes

::: gropt.diffusion.apply_recipe

::: gropt.diffusion.iter_recipes

::: gropt.diffusion.verify

::: gropt.diffusion.te_search

::: gropt.diffusion.sweep

::: gropt.diffusion.sweep_points

::: gropt.diffusion.to_raster

::: gropt.diffusion.continuation

::: gropt.diffusion.dt_schedule

::: gropt.diffusion.DiffParams
    options:
        summary: true
        show_signature_annotations: true

::: gropt.diffusion.SolverCfg
    options:
        summary: true
        show_signature_annotations: true

::: gropt.diffusion.SafeSource
    options:
        summary: true
        show_signature_annotations: true

::: gropt.diffusion.safe_ops

::: gropt.diffusion.build_gparams

::: gropt.diffusion.make_solver

::: gropt.diffusion_recipes.save_recipe

::: gropt.diffusion_recipes.load_recipe
