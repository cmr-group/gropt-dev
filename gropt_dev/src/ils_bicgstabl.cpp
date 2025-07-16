#include "spdlog/spdlog.h"

#include "ils_bicgstabl.hpp"

namespace Gropt {

ILS_BiCGstabl::ILS_BiCGstabl(GroptParams &_gparams, double _tol, double _sigma, int _n_iter, double _tik_lam)
    : IndirectLinearSolver(_gparams, _n_iter, _sigma, _tik_lam),
    tol(_tol)
{
    name = "BiCGStabl"; 

    int n = gparams->N * gparams->Naxis;
    rs.resize(ell + 1);
    us.resize(ell + 1);
    for (int i = 0; i <= ell; ++i) {
        rs[i].setZero(n);
        us[i].setZero(n);
    }

    Ax.setZero(n);
    b.setZero(n);

    tau.setZero(ell+1, ell+1);

    sigma.setZero(ell+1);
    gamma.setZero(ell+1);
    gammap.setZero(ell+1);
    gammapp.setZero(ell+1);
}


Eigen::VectorXd ILS_BiCGstabl::solve(Eigen::VectorXd &x_in)
{
    start_time = std::chrono::steady_clock::now();
    spdlog::trace("ILS_BiCGstabl::solve  start");

    x = x_in;

    get_rhs(x, b);
    get_lhs(x, Ax);

    r = (b - Ax);
    double rnorm = r.norm();
    double rnorm0 = rnorm;

    r_shadow = r;  // This could be random too.

    double omega = 1.0;
    double rho0 = 1.0;
    double alpha = 0.0;

    for (int i = 0; i <= ell; ++i) {
        us[i].setZero();
        rs[i].setZero();
    }
    rs[0] = r;

    int ii;
    for (ii = 0; ii < n_iter; ii++) {
        // First row from algorithm handled above, and redundant with last line (in printed algorithm)
        
        rho0 = -omega*rho0;

        for (int jj = 0; jj < ell; jj++) {
            double rho1 = r_shadow.dot(rs[jj]);
            double beta = alpha*rho1/rho0;
            rho0 = rho1;
            for (int ij = 0; ij <= jj; ij++) {
                us[ij] = rs[ij] - beta * us[ij];
            }
            get_lhs(us[jj], us[jj+1]);
            double gamma_cg = us[jj+1].dot(r_shadow);
            alpha = rho0/gamma_cg;
            for (int ij = 0; ij <= jj; ij++) {
                rs[ij] = rs[ij] - alpha * us[ij+1];
            }
            get_lhs(rs[jj], rs[jj+1]);
            x += alpha * us[0];
        }

        for (int jj = 1; jj <= ell; jj++) {
            for (int ij = 1; ij < jj; ij++) {
                tau(ij, jj) = rs[jj].dot(rs[ij])/sigma(ij);
                rs[jj] = rs[jj] - tau(ij, jj) * rs[ij];
            }
            sigma(jj) = rs[jj].dot(rs[jj]);
            gammap(jj) = rs[0].dot(rs[jj])/sigma(jj);
        }

        gamma(ell) = gammap(ell);
        omega = gamma(ell);

        for (int jj = ell-1; jj > 0; jj--) {
            double sum = 0.0;
            for (int ij = jj+1; ij <= ell; ij++) {
                sum += tau(jj, ij) * gamma(ij);
            }
            gamma(jj) = gammap(jj) - sum;
        }

        for (int jj = 1; jj < ell; jj++) {
            double sum = 0.0;
            for (int ij = jj+1; ij < ell; ij++) {
                sum += tau(jj, ij) * gamma(ij+1);
            }
            gammapp(jj) = gamma(jj+1) + sum;
        }

        x += gamma(1) * rs[0];
        rs[0] = rs[0] - gammap(ell)*rs[ell];
        us[0] = us[0] - gamma(ell)*us[ell];

        for (int jj = 1; jj < ell; jj++) {
            us[0] = us[0] - gamma(jj)*us[jj];
            x += gammapp(jj) * rs[jj];
            rs[0] = rs[0] - gammap(jj)*rs[jj];
        }
        rnorm = rs[0].norm();
        if ((rnorm <= tol * rnorm0))
        {
            spdlog::trace("ILS_CG::solve  break for (res <= tol)  ii = {:d}", ii);
            break;
        }
    }

    spdlog::info("ILS_BiCGstabl::solve  rnorm0 = {:e}   rnorm = {:e}   ii = {:d}", rnorm0, rs[0].norm(), ii);

    stop_time = std::chrono::steady_clock::now();
    elapsed_us = stop_time - start_time;

    hist_n_iter.push_back(ii+1);

    return x;
}


} // namespace Gropt