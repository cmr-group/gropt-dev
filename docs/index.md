---
hide:
  - navigation
  - toc
---

<div class="gropt-hero" markdown>

![GrOpt](img/gropt_logo.png){ .gropt-hero__logo }

<h1 class="gropt-hero__tagline">A toolbox for designing optimal MRI gradient waveforms</h1>

<p class="gropt-hero__sub">Rather than combinations of trapezoids, GrOpt uses constrained optimization to find the
time-optimal gradient waveform for a given set of hardware limits and encoding requirements. The solver is written
in C++, with Python bindings, and is designed to be compiled into a pulse sequence.</p>

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

## What you can design

Each card links to a worked example.

<div class="grid cards" markdown>

-   :material-sine-wave:{ .lg .middle } __Gradient moments__

    ---

    Null or set moments of any order, for flow compensation, phase encoding or rephasing, and find the
    shortest waveform that meets them.

    [:octicons-arrow-right-24: Example: flow encoding](examples/moments/flow_encoding.ipynb)

-   :material-pin-outline:{ .lg .middle } __Fixed gradient values__

    ---

    Hold any samples at set values: nonzero end points that join neighbouring gradients, or a slice-select
    plateau with its refocusing lobe designed around it.

    [:octicons-arrow-right-24: Example: fixed gradient values](examples/moments/nonzero_ends.ipynb)

-   :material-vector-curve:{ .lg .middle } __Waveform shape__

    ---

    Choose among the waveforms that meet the constraints: minimum gradient energy, or a total-variation
    penalty for blockier, trapezoid-like lobes.

    [:octicons-arrow-right-24: Example: regularization](examples/objective_functions/regularizations.ipynb)

-   :material-flash-outline:{ .lg .middle } __Nerve stimulation__

    ---

    Peripheral nerve and cardiac stimulation limits with the SAFE model, from your scanner's parameters.

    [:octicons-arrow-right-24: Example: PNS limits](examples/hardware/pns.ipynb)

-   :material-axis-arrow:{ .lg .middle } __Hardware limits__

    ---

    Gradient amplitude and slew rate across several axes, per axis or rotation invariant, and limits that change
    over time.

    [:octicons-arrow-right-24: Example: multiple axes](examples/hardware/multi_axis.ipynb)

-   :material-water-outline:{ .lg .middle } __Diffusion encoding__

    ---

    Maximize the b-value or find the minimum TE, with motion, eddy-current and concomitant-field compensation,
    from a short problem description.

    [:octicons-arrow-right-24: Example: diffusion encoding](examples/diffusion/diffusion.ipynb)

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

If GrOpt is useful in your work, please cite:

Loecher M, Middione MJ, Ennis DB. A gradient optimization toolbox for general purpose time-optimal MRI gradient
waveform design. *Magn Reson Med.* 2020; 84: 3234–3245.
[doi:10.1002/mrm.28384](https://doi.org/10.1002/mrm.28384)

??? quote "BibTeX"

    ```bibtex
    @article{loecher2020gropt,
      author  = {Loecher, Michael and Middione, Matthew J. and Ennis, Daniel B.},
      title   = {A gradient optimization toolbox for general purpose time-optimal {MRI} gradient waveform design},
      journal = {Magnetic Resonance in Medicine},
      volume  = {84},
      number  = {6},
      pages   = {3234--3245},
      year    = {2020},
      doi     = {10.1002/mrm.28384},
      url     = {https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.28384}
    }
    ```
