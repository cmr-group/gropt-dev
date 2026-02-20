from . import gropt_wrapper

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
