# Utilities

## Waveform Metrics and Plots

`gropt.plot` computes moments, b-values and the concomitant balance of a waveform, and plots waveforms.

::: gropt.plot.get_moments
    options:
        heading_level: 3

::: gropt.plot.get_bval
    options:
        heading_level: 3

::: gropt.plot.get_concomitant
    options:
        heading_level: 3

::: gropt.plot.plot_waves
    options:
        heading_level: 3

::: gropt.plot.plot_diff
    options:
        heading_level: 3

## Raster Conversion

A waveform solved on a coarse raster can be moved exactly onto the hardware raster, since the scanner plays the
linear interpolant between samples. For diffusion results, `gropt.diffusion.to_raster` does this in one call.

::: gropt.resample_waveform
    options:
        heading_level: 3

::: gropt.resample_inv_vec
    options:
        heading_level: 3

::: gropt.pwl_moment_weight
    options:
        heading_level: 3

::: gropt.low_freq_project
    options:
        heading_level: 3

## Logging

The C++ core logs through spdlog, at the warning level by default. In Jupyter its messages only reach the notebook
when forwarded with `setup_logging`.

::: gropt.utils.setup_logging
    options:
        heading_level: 3

::: gropt.set_log_level
    options:
        heading_level: 3
