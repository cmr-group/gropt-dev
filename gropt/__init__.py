from importlib.metadata import version

__version__ = version('gropt')

from . import debug_plotting
from . import readasc
from .gropt_wrapper import *
from .gropt_wrapper import __build_date__
from .readasc import get_random_safe_params
from .utils import demo, setup_logging

# `gropt.diff` was renamed to `gropt.diffusion`. The old name still works via the gropt/diff.py
# shim, which forwards to gropt.diffusion and emits a DeprecationWarning on use.
def __getattr__(name):
    if name == "diff":
        import importlib
        return importlib.import_module("gropt.diff")
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


set_log_level(3)  # warn
