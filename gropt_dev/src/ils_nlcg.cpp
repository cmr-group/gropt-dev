#include "spdlog/spdlog.h"

#include "ils_nlcg.hpp"

namespace Gropt {



ILS_NLCG::ILS_NLCG(GroptParams &_gparams)
    : IndirectLinearSolver(_gparams)
{
    name = "NLCG"; 

    b.setZero(gparams->N);
    Ax.setZero(gparams->N);
    x.setZero(gparams->N);
}


double ILS_NLCG::line_search()
{
    for(int i = -4; i < 4; i++) {
        double a = std::pow(10.0, i);
        
        Eigen::VectorXd xadx = x + a*dx;
        Ax.setZero();
        get_lhs(xadx, Ax);
    }
}

Eigen::VectorXd ILS_NLCG::solve(Eigen::VectorXd &x0)
{
    start_time = std::chrono::steady_clock::now();
    spdlog::trace("ILS_CG::solve  start");

    x = x0;

    b.setZero();
    get_rhs(x0, b);

    Ax.setZero();
    get_lhs(x, Ax);

    r = (Ax - b);
    double rnorm0 = r.norm();

    dx.setZero();
    get_lhs(r, dx);





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

    // tol0 = std::max(rel_tol*rnorm0/bnorm0, 1.0e-12);
    // if (gparams->iiter > 3) {
    //     tol = std::min(tol0, tol);
    // } else {
    //     tol = tol0;
    // }

    tol = 1e-12;

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

        if ((gamma <= tol * rnorm0) && (ii > 4))
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