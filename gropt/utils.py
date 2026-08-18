import atexit
import logging
from . import gropt_wrapper

# Maps spdlog levels (0–5) to Python logging levels.
# spdlog: 0=trace, 1=debug, 2=info, 3=warn, 4=err, 5=critical
_SPDLOG_TO_PYTHON = [
    logging.DEBUG,    # 0: trace
    logging.DEBUG,    # 1: debug
    logging.INFO,     # 2: info
    logging.WARNING,  # 3: warning
    logging.ERROR,    # 4: error
    logging.CRITICAL, # 5: critical
]


def _in_jupyter() -> bool:
    """True under a Jupyter/IPython kernel, where the raw C++ stdout is NOT shown in cells."""
    try:
        from IPython import get_ipython

        return get_ipython().__class__.__name__ == "ZMQInteractiveShell"
    except Exception:
        return False


def setup_logging(level: int = 2, to_python: bool | None = None) -> None:
    """Set the gropt C++ log level, optionally routing messages through Python's ``logging``.

    In a terminal the C++ logger already writes to stdout/stderr, so only the level matters and this just
    calls :func:`set_log_level` -- no Python callback, nothing to go wrong at exit. In Jupyter the raw C++
    stdout is not shown in cells, so messages are instead bridged through Python's ``logging`` (a
    ``StreamHandler``); that bridge is auto-enabled only when a Jupyter/IPython kernel is detected.

    Call once at the top of your script or notebook; repeat calls just update the level.

    Parameters
    ----------
    level : int
        spdlog level (lower = more verbose): 0=Trace, 1=Debug, 2=Info, 3=Warning, 4=Error, 5=Critical,
        6=Off. Default is 2 (Info).
    to_python : bool or None
        Route C++ logs through Python's ``logging`` module. ``None`` (default) auto-detects Jupyter.
        ``False`` keeps the native stdout/stderr sink -- the right choice in a terminal. ``True`` forces
        the bridge (needed to see logs in notebook cells).
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
            # The C++ default logger holds the Python callback above; release it while the interpreter is
            # still alive, or its destruction during C++ static teardown (after Py_Finalize) segfaults at exit.
            atexit.register(gropt_wrapper.clear_log_callback)

    gropt_wrapper.set_log_level(level)


def demo(plot=False):
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
