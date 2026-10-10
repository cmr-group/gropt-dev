::: gropt.GroptParams
    options:
        members: false
        heading_level: 2

## Helper Functions

These helper functions set all necessary parameters to define the problem.  They define some common usage patterns for optimization, such as fully optimized waveforms that start and end at 0 (`vec_init_simple`), or the optimization of diffusion waveforms (`diff_init`).

::: gropt.GroptParams.vec_init_simple

::: gropt.GroptParams.diff_init

::: gropt.GroptParams.diff_init_deadtime

::: gropt.GroptParams.diff_init_preencode

## Direct Setup

The individual components of the problem definition are:

| Name | Description |
| ---- | ----------- |
| `N` | The number of points of the waveform, per-axis |
| `Naxis` |  Number of axes |
| `dt` |  Gradient raster time [seconds] |
| `set_vals` | A vector that contains either gradient values or NaN, which describes if a given point in the waveform has a fixed value, or is optimizable (NaN) |
| `inv_vec` | A vector of either 1 or -1 describing which points of the waveform have been RF inverted |
| `X0` | An initial guess at the solution waveform |

Many of these will be automatically generated for the simplest options if omitted.

Vectors are flat and axis-major: length `Naxis * N`, with `X.reshape(Naxis, N)` giving one row per axis.

### `setvec_X0`

```python
setvec_X0(X0, set_others=True)
```

Set the initial waveform guess.

- **X0** (*np.ndarray*): initial guess [T/m], length `N` (sets `Naxis = 1`) or shape `(Naxis, N)`.
- **set_others** (*bool*, optional): if True, reset `inv_vec` to +1 and fix each axis's first and last points to
  their `X0` values (all others free). Use False to keep a layout from `diff_init`.

### `setvec_set_vals`

```python
setvec_set_vals(set_vals)
```

Set the fixed-value vector and rebuild the fixer mask. `X0` is set to the fixed values at fixed points; `inv_vec`
is not changed, so call `vec_init_simple()` or a `diff_init` variant first (see the
[fixed values example](../examples/moments/nonzero_ends.ipynb)).

- **set_vals** (*np.ndarray*): length `N` (sets `Naxis = 1`) or shape `(Naxis, N)`. NaN = free, finite value =
  fixed to that value.

::: gropt.GroptParams.getvec_set_vals

::: gropt.GroptParams.getvec_fixer

::: gropt.GroptParams.getvec_inv_vec

::: gropt.GroptParams.getvec_X0

::: gropt.GroptParams.set_ils_solver

## Operators

Each `add_*` call adds an operator. `solve()` prepares them automatically; these functions inspect them.

::: gropt.GroptParams.prepare

::: gropt.GroptParams.get_op_names

::: gropt.GroptParams.get_op_keys

::: gropt.GroptParams.print_op_details

::: gropt.GroptParams.reset_op_weights
