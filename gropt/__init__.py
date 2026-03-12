from importlib.metadata import version
__version__ = version("gropt")

from .gropt_wrapper import *
from .gropt_wrapper import __build_date__
from . import readasc
from .utils import demo, setup_logging

set_log_level(3)  # warn
