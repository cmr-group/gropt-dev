import gropt


def diff_solve(gparams, extra_iters=2000, ils_max_iter=30):
    solver = gropt.SolverGroptSDMM()
    solver.set_general_params(max_feval=200000, max_iter=20000, gamma_x=1.5, extra_iters=extra_iters)
    solver.set_ils_params(ils_max_iter=ils_max_iter, ils_tol=1e-12, ils_sigma=0.0001, ils_tik_lam=0.0001)
    result = solver.solve(gparams)
    return result
