#include <algorithm>

#include "spdlog/spdlog.h"

#include "ils_cg.hpp"

namespace Gropt {

ILS_CG::ILS_CG(GroptParams &_gparams, double _tol, int _min_iter, double _sigma, int _n_iter, double _tik_lam)
    : IndirectLinearSolver(_gparams, _n_iter, _sigma, _tik_lam), tol(_tol), min_iter(_min_iter) {
    name = "CG";

    int size = gparams->N * gparams->Naxis;
    b.setZero(size);
    Ax.setZero(size);
    Ap.setZero(size);
    r.setZero(size);
    p.setZero(size);
    x.setZero(size);
}

Eigen::VectorXd ILS_CG::solve(Eigen::VectorXd &x0) {
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
    // These were debugging variables that didn't prove particularly useful, could consider removing
    bool neg_curv = false;
    double min_curv = 0.0;
    double max_curv = 0.0;
    bool curv_init = false;

    x = x0;
    if (gparams->eq_proj.active) {
        gparams->eq_proj.project_affine(x); // enforce equality constraints exactly.
    }
    Eigen::VectorXd x_out = x;
    double r_min = std::numeric_limits<double>::max();

    b.setZero();
    get_rhs(x0, b);

    Ax.setZero();
    Ap.setZero();
    get_lhs(x, Ax);

    r = (b - Ax);
    if (gparams->eq_proj.active) {
        gparams->eq_proj.project_dir(r); // keep CG residual/directions in the equality null-space
    }
    rnorm0 = r.norm();
    bnorm0 = b.norm();
    // Warm-start-relative reduction (auto-tightens as ADMM converges) floored by a ||b||-relative
    // absolute tolerance (prevents an unreachable target when rnorm0 -> 0).
    double stop_thresh = std::max(tol * rnorm0, tol_abs_rel * bnorm0);

    p = r;
    gamma = r.dot(r);

    double gamma_new;
    int ii;
    for (ii = 0; ii < n_iter; ii++) {
        spdlog::trace("ILS_CG::solve  ii = {:d}  start", ii);

        Ap.setZero();
        get_lhs(p, Ap); // Ap = A*p
        if (gparams->eq_proj.active) {
            gparams->eq_proj.project_dir(Ap); // projected CG: effective operator is P M P
        }
        pAp = p.dot(Ap);
        
        if (pAp <= 0.0) {
            neg_curv = true; // non-positive curvature: indefinite/singular system
        }

        if (collect_debug) {  // Ignored if extra_debug is off
            double pp_curv = p.dot(p); // extra O(N) dot — only when debug logging is on
            if (pp_curv > 0.0) {
                double rq = pAp / pp_curv; // Rayleigh quotient = curvature of M along p
                if (!curv_init) {
                    min_curv = rq;
                    max_curv = rq;
                    curv_init = true;
                } else {
                    if (rq < min_curv) min_curv = rq;
                    if (rq > max_curv) max_curv = rq;
                }
            }
        }
        alpha = gamma / pAp;

        x += alpha * p;
        r -= alpha * Ap;

        gamma_new = r.dot(r);
        beta = gamma_new / gamma;
        gamma = gamma_new;

        p = beta * p + r;

        if ((std::sqrt(gamma) <= stop_thresh) && (ii > min_iter)) {
            spdlog::trace("ILS_CG::solve  break for (res <= tol)  ii = {:d}", ii);
            break;
        }
    }

    spdlog::trace("ILS_CG::solve  rnorm0 = {:e}   rnorm = {:e}  gamma = {:e}  ii = {:d}", rnorm0, r.norm(), gamma, ii);

    stop_time = std::chrono::steady_clock::now();
    elapsed_us = stop_time - start_time;

    hist_n_iter.push_back(ii + 1);
    hist_rnorm0.push_back(rnorm0);
    hist_rnorm.push_back(r.norm());
    hist_bnorm0.push_back(bnorm0);
    hist_neg_curv.push_back(neg_curv ? 1 : 0);
    if (collect_debug) {
        hist_min_curv.push_back(min_curv);
        hist_max_curv.push_back(max_curv);
    }

    return x;
}

} // namespace Gropt
