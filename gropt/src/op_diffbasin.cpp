#include <algorithm>
#include <cmath>

#include "spdlog/spdlog.h"

#include "op_diffbasin.hpp"

namespace Gropt {

Op_DiffBasin::Op_DiffBasin(const ProblemData &_pdata, double _window_time, double _eps_factor, double _gmax,
                           double _weight_mod, bool _same_sign)
    : Operator(_pdata) {
    name = "DiffBasin";
    window_time = _window_time;
    eps_factor = _eps_factor;
    gmax = _gmax;
    weight_mod = _weight_mod;
    same_sign = _same_sign;
}

void Op_DiffBasin::init() {
    spdlog::trace("Op_DiffBasin::init");

    Ax_size = 2; // [mean pre-180 gradient, mean post-180 gradient]
    const int Nfull = pdata->Naxis * pdata->N;
    const int w = std::max(1, (int)std::round(window_time / pdata->dt));

    a_pre.setZero(Nfull);
    a_post.setZero(Nfull);

    // Single 180 in the middle: located where inv_vec flips from + to -.
    int i_flip = 0;
    while (i_flip < Nfull && pdata->inv_vec(i_flip) >= 0) i_flip++;

    // Post window: first w FREE samples at/after the flip. Walking outward and taking only free
    // samples automatically skips the fixed 180 RF block.
    int n_post = 0;
    for (int i = i_flip; i < Nfull && n_post < w; i++) {
        if (pdata->fixer(i) > 0.0) {
            a_post(i) = 1.0;
            n_post++;
        }
    }
    if (n_post > 0) a_post /= (double)n_post; // -> mean over the collected free samples

    // Pre window: last w FREE samples before the flip.
    int n_pre = 0;
    for (int i = i_flip - 1; i >= 0 && n_pre < w; i--) {
        if (pdata->fixer(i) > 0.0) {
            a_pre(i) = 1.0;
            n_pre++;
        }
    }
    if (n_pre > 0) a_pre /= (double)n_pre;

    if (n_pre == 0 || n_post == 0) {
        spdlog::warn("Op_DiffBasin: no free samples in the {} window (inv_vec flip at index {}); that "
                     "side will be left unconstrained. Check the 180 location / window_time.",
                     (n_pre == 0 ? "pre-180" : "post-180"), i_flip);
    }

    eps = eps_factor * gmax;

    // Preconditioner only (the prox un/re-normalizes exactly, so this does not change the physical
    // threshold): largest row norm of the unit-mean functionals.
    double n2 = std::max(a_pre.squaredNorm(), a_post.squaredNorm());
    spec_norm2 = (n2 > 0.0) ? n2 : 1.0;
    spec_norm = std::sqrt(spec_norm2);

    target = 0.0;
    tol0 = 0.0;
    tol = 0.0;

    if (do_init_weights) {
        obj_weight = 1.0;
        obj_weight *= weight_mod;
    }

    Operator::init();

    spdlog::trace("Op_DiffBasin: w={} n_pre={} n_post={} eps={:.3e}", w, n_pre, n_post, eps);
}

void Op_DiffBasin::forward(Eigen::VectorXd &X, Eigen::VectorXd &out) {
    out(0) = a_pre.dot(X);  // mean gradient, pre-180 window
    out(1) = a_post.dot(X); // mean gradient, post-180 window
}

void Op_DiffBasin::transpose(Eigen::VectorXd &X, Eigen::VectorXd &out) {
    out = a_pre * X(0) + a_post * X(1); // adjoint of forward (X is the 2-vector in Ax-space)
}

void Op_DiffBasin::prox(Eigen::VectorXd &X) {
    spdlog::trace("Starting Op_DiffBasin::prox");

    if (do_equil) {
        X.array() /= eq_rows.array();
    }
    X.array() *= spec_norm; // -> physical mean-gradient units

    // One-sided (half-space) projection: pre >= +eps always; post <= -eps (opposite/flip, default) or
    // post >= +eps (same_sign). Only a side that actually has a window (norm > 0) is enforced.
    if (a_pre.squaredNorm() > 0.0 && X(0) < eps) X(0) = eps;
    if (a_post.squaredNorm() > 0.0) {
        if (same_sign) {
            if (X(1) < eps) X(1) = eps;     // post >= +eps  (same sign as pre)
        } else {
            if (X(1) > -eps) X(1) = -eps;   // post <= -eps  (opposite / flip)
        }
    }

    if (do_equil) {
        X.array() *= eq_rows.array();
    }
    X.array() /= spec_norm;

    spdlog::trace("Finished Op_DiffBasin::prox");
}

void Op_DiffBasin::check(Eigen::VectorXd &X) {
    int is_feas = 1;

    if (do_equil) {
        X.array() /= eq_rows.array();
    }
    X.array() *= spec_norm;

    if (a_pre.squaredNorm() > 0.0 && X(0) < eps) is_feas = 0;
    if (a_post.squaredNorm() > 0.0) {
        if (same_sign) {
            if (X(1) < eps) is_feas = 0;
        } else {
            if (X(1) > -eps) is_feas = 0;
        }
    }

    if (do_equil) {
        X.array() *= eq_rows.array();
    }
    X.array() /= spec_norm;

    hist_feas.push_back(is_feas);
}

} // namespace Gropt
