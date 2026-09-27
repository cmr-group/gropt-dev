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

`converged` alone does not prove a waveform is within its limits: after a numerical blow-up the solver can return
an earlier iterate. `verify(result, cfg)` recomputes every limit from the waveform (PNS/CNS on the 10 us hardware
raster) and returns the worst ratio per limit and `feasible`.

## API

::: gropt.diffusion.solve
    options:
        heading_level: 3

::: gropt.diffusion.available_recipes
    options:
        heading_level: 3

::: gropt.diffusion.apply_recipe
    options:
        heading_level: 3

::: gropt.diffusion.iter_recipes
    options:
        heading_level: 3

::: gropt.diffusion.verify
    options:
        heading_level: 3

::: gropt.diffusion.te_search
    options:
        heading_level: 3

::: gropt.diffusion.DiffParams
    options:
        heading_level: 3

::: gropt.diffusion.SolverCfg
    options:
        heading_level: 3

::: gropt.diffusion_recipes.save_recipe
    options:
        heading_level: 3

::: gropt.diffusion_recipes.load_recipe
    options:
        heading_level: 3
