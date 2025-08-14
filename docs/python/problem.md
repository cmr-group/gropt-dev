## Helper Functions

These helper functions set all necessary parameterss to define the function.

::: gropt.GroptParams
    options:
        members:
            - vec_init_simple
            - diff_init
        heading_level: 3

## Direct Setup

The individual components of the problem definition are:

| Name | Description |
| ---- | ----------- |
| `N` | The number of points of the wave form, per-axis |
| `Naxis` |  Number of axis |
| `dt` |  Gradient raster time [seconds] |
| `set_vec` | A vector that contains either gradient values or NaN, which describes if a given point in the waveform has a fixed value, or is optimizable (NaN) |
| `inv_vec` | A vector of either 1 or -1 describing which points of the waveform have been RF inverted |
| `X0` | An initial guess at the solution waveform |

Many of these wll be automatically generated for the simplest options if omitted.

::: gropt.GroptParams
    options:
        members:
            - setvec_X0
        heading_level: 3