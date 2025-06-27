from . import gropt_wrapper
import matplotlib.pyplot as plt

def demo():
    print('Starting demo...', flush = True)
    
    gparams = gropt_wrapper.GroptParams()
    gparams.N = 102
    gparams.Naxis = 1
    gparams.dt = 10e-6
    gparams.vec_init_simple()

    gparams.add_gmax(.08)
    gparams.add_smax(200)
    gparams.add_moment(0,2.0)
    gparams.add_moment(1,0.0)
    gparams.add_moment(2,0.0)
    
    print('Starting solve...', flush = True)

    gparams.init()
    gparams.solve()

    print('Finished solve...', flush = True)

    out = gparams.get_out()
    
    print(f'{out.shape = }')

    # plt.figure()
    # plt.plot(out)
    # plt.show()
    
    print('Done!', flush = True)
