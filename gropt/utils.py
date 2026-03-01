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


def setup_logging(level: int = 2) -> None:
    """Route gropt C++ log messages through Python's logging module.

    Call once at the top of your script or notebook. Subsequent calls
    update the log level without adding duplicate handlers.

    Parameters
    ----------
    level : int
        spdlog level (lower = more verbose): 0=Trace, 1=Debug, 2=Info,
        3=Warning, 4=Error, 5=Critical, 6=Off. Default is 2 (Info).
    """
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
