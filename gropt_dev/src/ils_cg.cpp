#include "spdlog/spdlog.h"

#include "ils_cg.hpp"

namespace Gropt {

ILS_CG::ILS_CG(GroptParams &_gparams, double _tol, int _min_iter, double _sigma, int _n_iter, double _tik_lam)
    : IndirectLinearSolver(_gparams, _n_iter, _sigma, _tik_lam), 
    tol(_tol), 
    min_iter(_min_iter)
{
    name = "CG"; 

    b.setZero(gparams->N);
    Ax.setZero(gparams->N);
    Ap.setZero(gparams->N);
    r.setZero(gparams->N);
    p.setZero(gparams->N);
    x.setZero(gparams->N);
}

Eigen::VectorXd ILS_CG::solve(Eigen::VectorXd &x0)
{
    start_time = std::chrono::steady_clock::now();
    spdlog::trace("ILS_CG::solve  start");

    double rnorm0;
    double bnorm0;
    double tol0;
    double res;

    double alpha; 
    double beta;
    double gamma;

    double pAp;

    x = x0;

    b.setZero();
    get_rhs(x0, b);

    Ax.setZero();
    Ap.setZero();
    get_lhs(x, Ax);

    r = (b - Ax);
    rnorm0 = r.norm();
    bnorm0 = b.norm();

    p = r;
    gamma = r.dot(r);

    double gamma_new;
    int ii;
    for (ii = 0; ii < n_iter; ii++) {
        spdlog::trace("ILS_CG::solve  ii = {:d}  start", ii);

        Ap.setZero();
        get_lhs(p, Ap);  // Ap = A*p
        pAp = p.dot(Ap);
        alpha = gamma / pAp;

        x += alpha * p;
        r -= alpha * Ap;

        gamma_new = r.dot(r);
        beta = gamma_new / gamma;
        gamma = gamma_new;

        p = beta * p + r;

        if ((gamma <= tol * rnorm0) && (ii > min_iter))
        {
            spdlog::trace("ILS_CG::solve  break for (res <= tol)  ii = {:d}", ii);
            break;
        }

        
        
    }

    spdlog::trace("ILS_CG::solve  rnorm0 = {:e}   rnorm = {:e}   ii = {:d}", rnorm0, r.norm(), ii);

    stop_time = std::chrono::steady_clock::now();
    elapsed_us = stop_time - start_time;

    hist_n_iter.push_back(ii+1);

    return x;
}


} // namespace Gropt