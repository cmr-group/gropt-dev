## Helper Functions

These helper functions set all necessary parameters to define the problem.  They define some common usage patterns for optimization, such as fully optimized waveforms that start and end at 0 (`vec_init_simple`), or the optimization of diffusion waveforms (`diff_init`).

::: gropt.GroptParams
    options:
        members:
            - vec_init_simple
            - diff_init
            - diff_init_deadtime
            - diff_init_preencode
        heading_level: 3

## Direct Setup

The individual components of the problem definition are:

| Name | Description |
| ---- | ----------- |
| `N` | The number of points of the waveform, per-axis |
| `Naxis` |  Number of axes |
| `dt` |  Gradient raster time [seconds] |
| `set_vec` | A vector that contains either gradient values or NaN, which describes if a given point in the waveform has a fixed value, or is optimizable (NaN) |
| `inv_vec` | A vector of either 1 or -1 describing which points of the waveform have been RF inverted |
| `X0` | An initial guess at the solution waveform |

Many of these will be automatically generated for the simplest options if omitted.

::: gropt.GroptParams
    options:
        members:
            - setvec_X0
            - set_ils_solver
        heading_level: 3
