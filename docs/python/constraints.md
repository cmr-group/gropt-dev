## Gradient Amplitude

::: gropt.GroptParams
    options:
        members:
            - add_gmax
        heading_level: 3
        show_root_heading: false

## Slew Rate

::: gropt.GroptParams
    options:
        members:
            - add_smax
        heading_level: 3
        show_root_heading: false

## Moments

::: gropt.GroptParams
    options:
        members:
            - add_moment
        heading_level: 3
        show_root_heading: false

## b-value

::: gropt.GroptParams
    options:
        members:
            - add_bvalue
        heading_level: 3
        show_root_heading: false

## PNS (SAFE-model)

PNS parameters can be loaded from an .asc file using the `gropt.readasc.asc_to_safe` function, which will create a python dictionary that can be used as input to the SAFE constraints.  Alternatively, no parameters can be entered to use some default demonstration values.  Finally, to get a SAFE predcition, the `gropt.get_SAFE` function can be used.

::: gropt.readasc.asc_to_safe
    options:
        heading_level: 3

::: gropt.GroptParams
    options:
        members:
            - add_SAFE
            - add_SAFE_vec
        heading_level: 3
        show_root_heading: false

::: gropt.get_SAFE
    options:
        heading_level: 3

