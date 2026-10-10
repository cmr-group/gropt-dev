---
hide:
  - navigation
  - toc
---

<div class="gropt-hero" markdown>

![GrOpt](img/gropt_logo.png){ .gropt-hero__logo }

<h1 class="gropt-hero__tagline">Time-optimal MRI gradient waveforms, designed by constrained optimization</h1>

<p class="gropt-hero__sub">Describe what a gradient waveform has to do (moments, b-value, timing, hardware and nerve
stimulation limits) and GrOpt finds the shortest waveform that does it, from a C++ core built to run inside a pulse
sequence.</p>

[Get started](examples/simple_example.ipynb){ .md-button .md-button--primary }
[Examples](examples/moments/flow_encoding.ipynb){ .md-button }
[GitHub](https://github.com/cmr-group/gropt-dev){ .md-button }

`pip install --pre gropt`

</div>

## Same b-value, shorter echo time

<figure class="gropt-figure" markdown>
![Velocity-compensated diffusion encoding for b = 1000 s/mm²: a conventional symmetric design needs TE 86.4 ms, GrOpt 72.0 ms](img/hero_light.png#only-light)
![Velocity-compensated diffusion encoding for b = 1000 s/mm²: a conventional symmetric design needs TE 86.4 ms, GrOpt 72.0 ms](img/hero_dark.png#only-dark)
<figcaption markdown>Velocity-compensated diffusion encoding (M0 = M1 = 0) for b = 1000 s/mm² with 80 mT/m and
100 T/m/s and a 16 ms readout. The conventional design is symmetric about the 180; GrOpt uses the time on both
sides. See [minimum TE](examples/diffusion/diffusion_te.ipynb).</figcaption>
</figure>

## What it does

<div class="grid cards" markdown>

-   :material-speedometer:{ .lg .middle } __Hardware limits__

    ---

    Gradient amplitude and slew rate, per axis or rotation invariant, constant or varying over time.

    [:octicons-arrow-right-24: Constraints](python/constraints.md)

-   :material-sine-wave:{ .lg .middle } __Moments and fixed samples__

    ---

    Null or target gradient moments of any order over any window, and fix any samples: nonzero ends, RF
    gaps, plateaus.

    [:octicons-arrow-right-24: Flow encoding](examples/moments/flow_encoding.ipynb)

-   :material-flash-outline:{ .lg .middle } __Nerve stimulation__

    ---

    Peripheral nerve and cardiac stimulation limits with the SAFE model, from your scanner's `.asc` file.

    [:octicons-arrow-right-24: PNS limits](examples/hardware/pns.ipynb)

-   :material-water-outline:{ .lg .middle } __Diffusion__

    ---

    Maximize the b-value or find the minimum TE, with motion, eddy-current and concomitant-field compensation.

    [:octicons-arrow-right-24: Diffusion encoding](examples/diffusion/diffusion.ipynb)

-   :material-timer-outline:{ .lg .middle } __Fast__

    ---

    Typical waveforms solve in milliseconds, and warm starts make searches over durations cheap.

    [:octicons-arrow-right-24: Solver control](examples/solvers/solver_control.ipynb)

-   :material-language-python:{ .lg .middle } __Python and C++__

    ---

    A Python package with notebooks for every feature, on the same C++ solver a sequence can embed.

    [:octicons-arrow-right-24: API reference](python/gropt_py.md)

</div>

## Quick start

A flow-compensated phase encode: M0 = 14 mT·ms/m and M1 = 0 in 1.12 ms, within 80 mT/m and 200 T/m/s.

<div class="grid" markdown>

```python
import gropt

gparams = gropt.GroptParams()
gparams.dt = 10e-6           # raster time [s]
gparams.N = 112              # 1.12 ms
gparams.vec_init_simple()    # free, starting and ending at 0

gparams.add_gmax(0.08)       # [T/m]
gparams.add_smax(200)        # [T/m/s]
gparams.add_moment(0, 14.0)  # M0 [mT*ms/m]
gparams.add_moment(1, 0.0)   # M1

result = gropt.solve(gparams)
waveform = result.X          # [T/m]
```

![The resulting flow-compensated phase-encode waveform](img/quickstart_light.png#only-light)
![The resulting flow-compensated phase-encode waveform](img/quickstart_dark.png#only-dark)

</div>

[:octicons-arrow-right-24: The quick start as a notebook](examples/simple_example.ipynb), and the
[Python API reference](python/gropt_py.md).

## Citing GrOpt

<!-- TODO(citations): the references below are placeholders written from memory; replace them with the exact
     citations (and add volume, pages and DOI) before publishing. -->

If GrOpt is useful in your work, please cite:

- M. Loecher, M. J. Middione, D. B. Ennis. *A gradient optimization toolbox for general purpose time-optimal MRI
  gradient waveform design.* Magnetic Resonance in Medicine, 2020.
- E. Aliotta, H. H. Wu, D. B. Ennis. *Convex optimized diffusion encoding (CODE) gradient waveforms for minimum
  echo time and bulk motion-compensated diffusion-weighted MRI.* Magnetic Resonance in Medicine, 2017.

??? quote "BibTeX"

    ```bibtex
    @article{loecher2020gropt,
      author  = {Loecher, Michael and Middione, Matthew J. and Ennis, Daniel B.},
      title   = {A gradient optimization toolbox for general purpose time-optimal {MRI} gradient waveform design},
      journal = {Magnetic Resonance in Medicine},
      year    = {2020}
    }

    @article{aliotta2017code,
      author  = {Aliotta, Eric and Wu, Holden H. and Ennis, Daniel B.},
      title   = {Convex optimized diffusion encoding ({CODE}) gradient waveforms for minimum echo time and bulk
                 motion-compensated diffusion-weighted {MRI}},
      journal = {Magnetic Resonance in Medicine},
      year    = {2017}
    }
    ```
