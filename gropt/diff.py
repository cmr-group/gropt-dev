"""Deprecated alias for :mod:`gropt.diffusion`.

The module was renamed ``diff`` -> ``diffusion``. Anything accessed through
``gropt.diff`` is forwarded to ``gropt.diffusion`` and triggers a
DeprecationWarning the first time, so existing code keeps working while you
migrate. This shim will be removed in a future release.
"""

import warnings

from . import diffusion as _diffusion

_warned = False


def __getattr__(name):
    # PEP 562 module-level __getattr__: fires for every name pulled through
    # gropt.diff (incl. `from gropt.diff import X` and `gropt.diff.X`), so it
    # forwards public AND private names and warns once.
    global _warned
    if not _warned:
        warnings.warn(
            "gropt.diff has been renamed to gropt.diffusion; update your imports "
            "(the gropt.diff alias will be removed in a future release).",
            DeprecationWarning,
            stacklevel=2,
        )
        _warned = True
    return getattr(_diffusion, name)


def __dir__():
    return dir(_diffusion)
