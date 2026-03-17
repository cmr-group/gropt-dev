from importlib.metadata import version

__version__ = version('gropt')

from . import readasc
from .gropt_wrapper import *
from .gropt_wrapper import __build_date__
from .readasc import get_random_safe_params
from .utils import demo, setup_logging

set_log_level(3)  # warn
