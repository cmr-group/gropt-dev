#include "spdlog/spdlog.h"

#include "ils_nlcg.hpp"

namespace Gropt {

ILS_NLCG::ILS_NLCG(GroptParams &_gparams, double _sigma, int _n_iter, double _tik_lam)
    : IndirectLinearSolver(_gparams, _n_iter, _sigma, _tik_lam)
{
    name = "NLCG"; 

    b.setZero(gparams->N);
    Ax.setZero(gparams->N);
    x0.setZero(gparams->N);
}


Eigen::VectorXd ILS_NLCG::solve(Eigen::VectorXd &x_in)
{
    start_time = std::chrono::steady_clock::now();
    spdlog::trace("ILS_NLCG::solve  start");

    x0 = x_in;

    get_rhs(x0, b);

    get_lhs(x0, Ax);

    r = (Ax - b);
    double rnorm0 = r.squaredNorm();
    double rnorm_x0 = rnorm0;

    Eigen::VectorXd g0;
    Eigen::VectorXd g1;

    g0.setZero(r.size());
    g1.setZero(r.size());

    get_lhs(r, g0);
    g0 *= 2;

    Eigen::VectorXd d0 = -g0;
    Eigen::VectorXd d1;

    Eigen::VectorXd xad;
    Eigen::VectorXd y;
    Eigen::VectorXd ymid;

    double alpha0 = 2.0;
    double alpha;
    double rnorm_xad;
    double beta;
    
    int ii;
    for (ii = 0; ii < n_iter; ii++) {
        alpha = alpha0;

        xad = x0 + alpha*d0;
        get_lhs(xad, Ax);
        r = (Ax - b);
        rnorm_xad = r.squaredNorm();
        
        while (rnorm_xad > rnorm_x0 + eta*alpha*g0.dot(d0)) {
            alpha *= theta;
            xad = x0 + alpha*d0;
            get_lhs(xad, Ax);
            r = (Ax - b);
            rnorm_xad = r.squaredNorm();
            if (alpha < 1e-32) {
                spdlog::warn("ILS_NLCG::solve  alpha < 1e-32, stopping line search");
                break;
            } 
        }

        x1 = xad;
        get_lhs(r, g1);
        g1 *= 2;

        y = g1 - g0;

        // ymid = y - 2*d0*y.squaredNorm()/d0.dot(y);
        // beta = 1.0/(d0.dot(y)) * ymid.dot(g1);

        beta = g1.dot(y) / g0.squaredNorm();
        if (beta < 0.0) {
            beta = 0.0;
        }

        // beta = g1.squaredNorm() / g0.squaredNorm();

        d1 = -g1 + beta*d0;

        if (g1.dot(d1) >= 0.0) {
            d1 = -g1;
        }

        x0 = x1;
        g0 = g1;
        d0 = d1;
        rnorm_x0 = rnorm_xad;
        alpha0 = 4 * alpha;
    }

    spdlog::info("ILS_NLCG::solve  rnorm0 = {:e}   rnorm = {:e}   alpha0 = {:e}   ii = {:d}", sqrt(rnorm0), r.norm(), alpha0, ii);

    stop_time = std::chrono::steady_clock::now();
    elapsed_us = stop_time - start_time;

    hist_n_iter.push_back(ii+1);

    return x1;
}


} // namespace Gropt