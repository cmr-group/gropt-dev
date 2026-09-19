"""Logging setup and a small demo solve."""

import atexit
import logging
from . import gropt_wrapper

# spdlog level (index) -> Python logging level
_SPDLOG_TO_PYTHON = [
    logging.DEBUG,    # 0: trace
    logging.DEBUG,    # 1: debug
    logging.INFO,     # 2: info
    logging.WARNING,  # 3: warning
    logging.ERROR,    # 4: error
    logging.CRITICAL, # 5: critical
]


def _in_jupyter() -> bool:
    """Return True under a Jupyter kernel (ZMQInteractiveShell), where C++ stdout is not shown in cells."""
    try:
        from IPython import get_ipython

        return get_ipython().__class__.__name__ == "ZMQInteractiveShell"
    except Exception:
        return False


def setup_logging(level: int = 2, to_python: bool | None = None) -> None:
    """Set the gropt C++ log level, optionally forwarding messages to Python ``logging``.

    In Jupyter, C++ log messages are forwarded to Python logging (auto-detected); otherwise only the C++
    log level is set. Forwarded messages go to the ``gropt`` logger, which gets its own ``StreamHandler``
    and does not propagate. Call once at the top of a script or notebook; repeat calls only update the level.

    Parameters
    ----------
    level : int
        spdlog level (lower = more verbose): 0=Trace, 1=Debug, 2=Info, 3=Warning, 4=Error, 5=Critical,
        6=Off. Default is 2 (Info).
    to_python : bool or None
        Forward C++ logs to Python ``logging``. ``None`` (default) enables this only in Jupyter; ``False``
        keeps the native stdout/stderr sink; ``True`` forces forwarding.
    """
    if to_python is None:
        to_python = _in_jupyter()

    if to_python:
        logger = logging.getLogger("gropt")
        logger.setLevel(logging.DEBUG)
        logger.propagate = False

        if not logger.handlers:
            handler = logging.StreamHandler()
            handler.setFormatter(logging.Formatter("| {levelname:>8} |  {message}", style="{"))
            logger.addHandler(handler)
            gropt_wrapper.set_log_callback(
                lambda lvl, msg: logger.log(
                    _SPDLOG_TO_PYTHON[lvl] if lvl < len(_SPDLOG_TO_PYTHON) else logging.DEBUG,
                    msg,
                )
            )
            # release the Python callback before interpreter shutdown, or exit segfaults
            atexit.register(gropt_wrapper.clear_log_callback)

    gropt_wrapper.set_log_level(level)


def demo(plot=False):
    """Solve a small single-axis problem (gmax, smax, M0 = 2, M1 = M2 = 0) and print a summary.

    Parameters
    ----------
    plot : bool
        If True, plot the resulting waveform with matplotlib.
    """
    print('Starting demo...', flush=True)

    gparams = gropt_wrapper.GroptParams()
    gparams.N = 102
    gparams.Naxis = 1
    gparams.dt = 10e-6
    gparams.vec_init_simple()

    gparams.add_gmax(.08)
    gparams.add_smax(200)
    gparams.add_moment(0, 2.0)
    gparams.add_moment(1, 0.0)
    gparams.add_moment(2, 0.0)

    print('Starting solve...', flush=True)

    result = gropt_wrapper.solve(gparams)

    print('Finished solve...', flush=True)
    print(f'{result.converged = }', flush=True)
    print(f'{result.X.shape = }', flush=True)

    if plot:
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(result.X)
        plt.show()

    print('Done!', flush=True)
