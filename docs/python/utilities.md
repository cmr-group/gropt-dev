# Utilities

## Waveform Metrics and Plots

`gropt.plot` computes moments, b-values and the concomitant balance of a waveform, and plots waveforms.

::: gropt.plot.get_moments

::: gropt.plot.get_bval

::: gropt.plot.get_concomitant

::: gropt.plot.plot_waves

::: gropt.plot.plot_diff

## Raster Conversion

A waveform solved on a coarse raster can be moved exactly onto the hardware raster, since the scanner plays the
linear interpolant between samples. For diffusion results, `gropt.diffusion.to_raster` does this in one call.

::: gropt.resample_waveform

::: gropt.resample_inv_vec

::: gropt.pwl_moment_weight

::: gropt.low_freq_project

## Logging

The C++ core logs through spdlog, at the warning level by default. In Jupyter its messages only reach the notebook
when forwarded with `setup_logging`.

::: gropt.utils.setup_logging

::: gropt.set_log_level
