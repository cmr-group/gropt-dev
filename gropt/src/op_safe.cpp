#include <stdexcept>

#include "spdlog/spdlog.h"

#include "op_safe.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace Gropt {

Op_SAFE::Op_SAFE(const ProblemData &_pdata, double _stim_thresh, double _weight_mod) : Operator(_pdata) {
    name = "SAFE";
    stim_thresh = _stim_thresh;
    weight_mod = _weight_mod;
}

Op_SAFE::Op_SAFE(const ProblemData &_pdata, const Eigen::VectorXd &_stim_thresh_vec, double _weight_mod)
    : Operator(_pdata) {
    name = "SAFE";
    stim_thresh = 1.0;
    weight_mod = _weight_mod;
    stim_thresh_vec = _stim_thresh_vec;
}

void Op_SAFE::init() {
    spdlog::trace("Op_SAFE::init  N = {}", pdata->N);

    if (!rot_variant) {
        throw std::invalid_argument("Op_SAFE: rot_variant=false (rotationally invariant) is not implemented yet");
    }

    safe_params.calc_alphas(pdata->dt);

    if ((safe_params.a3[0] == 0) && (safe_params.a3[1] == 0) && (safe_params.a3[2] == 0)) {
        n_terms = 2;
        spdlog::trace("Op_SAFE::init  n_terms = {}", n_terms);
    }

    // target/tol0 are only logged; prox and check use stim_thresh_vec
    target = 0;
    tol0 = stim_thresh;
    tol = (1.0 - cushion) * tol0;

    if (stim_thresh_vec.size() != pdata->Naxis * pdata->N) {
        if (stim_thresh_vec.size() != 0) {
            spdlog::warn("Op_SAFE::init  stim_thresh_vec size does not match Naxis * N, resizing to match");
        }
        stim_thresh_vec.resize(pdata->Naxis * pdata->N);
        for (int i = 0; i < stim_thresh_vec.size(); i++) {
            stim_thresh_vec(i) = stim_thresh;
        }
    }

    Ax_size = n_terms * pdata->Naxis * pdata->N;

    signs1.setZero(pdata->Naxis * pdata->N);
    signs2.setZero(pdata->Naxis * pdata->N);
    signs3.setZero(pdata->Naxis * pdata->N);
    stim1.setZero(pdata->Naxis * pdata->N);
    stim2.setZero(pdata->Naxis * pdata->N);
    stim3.setZero(pdata->Naxis * pdata->N);
    slew_temp.setZero(pdata->Naxis * pdata->N);

    Operator::init();

    // No closed form; estimate ||A|| numerically
    spec_norm = estimate_self_spec_norm(30);
    spec_norm2 = spec_norm * spec_norm;

    spdlog::trace("Op_SAFE::init  spec_norm = {:.4f}", spec_norm);
    spdlog::trace("Op_SAFE::init  Done!");
}

void Op_SAFE::compute_slew(const Eigen::VectorXd &X) {
    for (int j = 0; j < Naxis; j++) {
        slew_temp(j * N) = X(j * N) / dt;
        for (int i = 1; i < N; i++) slew_temp(j * N + i) = (X(j * N + i) - X(j * N + i - 1)) / dt;
    }
}

void Op_SAFE::lowpass(Eigen::VectorXd &v, const std::vector<double> &alpha) const {
    for (int j = 0; j < Naxis; j++) {
        v(j * N) = alpha[j] * v(j * N);
        for (int i = 1; i < N; i++) v(j * N + i) = alpha[j] * v(j * N + i) + (1 - alpha[j]) * v(j * N + i - 1);
    }
}

void Op_SAFE::lowpass_T(Eigen::VectorXd &v, const std::vector<double> &alpha) const {
    for (int j = 0; j < Naxis; j++) {
        v(j * N + N - 1) = alpha[j] * v(j * N + N - 1);
        for (int i = N - 2; i >= 0; i--) v(j * N + i) = alpha[j] * v(j * N + i) + (1 - alpha[j]) * v(j * N + i + 1);
    }
}

void Op_SAFE::take_abs(Eigen::VectorXd &v, Eigen::VectorXd &signs) const {
    for (int i = 0; i < v.size(); i++) {
        const double x = v(i);
        if (freeze_signs) {
            v(i) = signs(i) * x; // frozen linearization: held sign, no recapture
        } else if (safe_eps > 0.0) {
            const double sa = sqrt(x * x + safe_eps * safe_eps);
            signs(i) = x / sa; // smooth sign in [-1,1], continuous through 0 (= d/dx of softabs)
            v(i) = sa;
        } else {
            signs(i) = (x < 0.0) ? -1.0 : 1.0;
            v(i) = (x < 0.0) ? -x : x;
        }
    }
}

void Op_SAFE::diff_T(Eigen::VectorXd &out) const {
    for (int j = 0; j < Naxis; j++) {
        for (int i = 0; i < N - 1; i++) out(j * N + i) = (out(j * N + i) - out(j * N + i + 1)) / dt;
        out(j * N + N - 1) = out(j * N + N - 1) / dt;
    }
}

void Op_SAFE::forward(Eigen::VectorXd &X, Eigen::VectorXd &out) {
    compute_slew(X);

    stim1 = slew_temp;                                   // term 1: |LP1(slew)|
    lowpass(stim1, safe_params.alpha1);
    if (!signed_terms13) take_abs(stim1, signs1);

    stim2 = slew_temp;                                   // term 2: LP2(|slew|), abs INSIDE the filter
    take_abs(stim2, signs2);
    lowpass(stim2, safe_params.alpha2);

    if (n_terms == 3) {                                  // term 3: |LP3(slew)|
        stim3 = slew_temp;
        lowpass(stim3, safe_params.alpha3);
        if (!signed_terms13) take_abs(stim3, signs3);
    }

    // Per axis j: out = [a1 stim1; a2 stim2; a3 stim3] * g_scale / stim_limit (n_terms blocks of N)
    for (int j = 0; j < Naxis; j++) {
        for (int i = 0; i < N; i++) {
            out(j * n_terms * N + i) =
                safe_params.a1[j] * stim1(j * N + i) / safe_params.stim_limit[j] * safe_params.g_scale[j];
            out(j * n_terms * N + i + N) =
                safe_params.a2[j] * stim2(j * N + i) / safe_params.stim_limit[j] * safe_params.g_scale[j];
            if (n_terms == 3) {
                out(j * n_terms * N + i + 2 * N) =
                    safe_params.a3[j] * stim3(j * N + i) / safe_params.stim_limit[j] * safe_params.g_scale[j];
            }
        }
    }
}

void Op_SAFE::transpose(Eigen::VectorXd &X, Eigen::VectorXd &out) {

    // Split and scale by a, stim_limit, g_scale
    for (int j = 0; j < Naxis; j++) {
        for (int i = 0; i < N; i++) {
            stim1(j * N + i) =
                safe_params.a1[j] * X(j * n_terms * N + i) / safe_params.stim_limit[j] * safe_params.g_scale[j];
            stim2(j * N + i) =
                safe_params.a2[j] * X(j * n_terms * N + i + N) / safe_params.stim_limit[j] * safe_params.g_scale[j];
            if (n_terms == 3) {
                stim3(j * N + i) = safe_params.a3[j] * X(j * n_terms * N + i + 2 * N) / safe_params.stim_limit[j] *
                                   safe_params.g_scale[j];
            }
        }
    }

    // term 1 applies its sign before the filter adjoint; term 2 after, mirroring the forward order
    if (!signed_terms13) stim1.array() *= signs1.array();
    lowpass_T(stim1, safe_params.alpha1);

    lowpass_T(stim2, safe_params.alpha2);
    stim2.array() *= signs2.array();

    if (n_terms == 3) {
        if (!signed_terms13) stim3.array() *= signs3.array();
        lowpass_T(stim3, safe_params.alpha3);
    }

    for (int j = 0; j < Naxis; j++) {
        for (int i = 0; i < N; i++) {
            if (n_terms == 3) {
                out(j * N + i) = stim1(j * N + i) + stim2(j * N + i) + stim3(j * N + i);
            } else if (n_terms == 2) {
                out(j * N + i) = stim1(j * N + i) + stim2(j * N + i);
            }
        }
    }

    diff_T(out); // out = D^T(stim1 + stim2 + stim3)
}

double Op_SAFE::axis_stim(const Eigen::VectorXd &X, int j, int i) const {
    double v = X(j * n_terms * N + i) + X(j * n_terms * N + i + N);
    if (n_terms == 3) v += X(j * n_terms * N + i + 2 * N);
    return v;
}

void Op_SAFE::prox_signed(Eigen::VectorXd &X) {
    if (do_equil) {
        X.array() /= eq_rows.array();
    }
    X.array() *= spec_norm;

    // With terms 1 and 3 emitted signed, the set is
    //     { y : sum_j ( (|y_j1| + |y_j2| + |y_j3|) / th_j )^2 <= upper^2 },
    // symmetric in every coordinate, so the projection keeps each sign and only shrinks magnitudes: work
    // with b = |y| and put the signs back at the end.
    //
    // With a multiplier lam >= 0 the KKT conditions shrink every term of an axis by the SAME amount
    // delta_j (the equal share of the metric projection), floored at zero:
    //     w_jk = max(b_jk - delta_j, 0),   delta_j = lam * s_j / th_j^2,   s_j = sum_k w_jk.
    // Over the terms still above the floor (the active set A_j) that closes:
    //     delta_j = lam B_j / (th_j^2 + lam n_j),   s_j = B_j th_j^2 / (th_j^2 + lam n_j),
    // with B_j the sum and n_j the count of the active terms. Newton on lam, re-resolving A_j each step;
    // with all terms active this reduces to Op_SAFE::prox's formula exactly.
    const double upper = 1.0 - cushion;
    const int nt = n_terms;
    std::vector<double> b(Naxis * nt), th(Naxis), B(Naxis), sgn(Naxis * nt);
    std::vector<int> nact(Naxis);
    std::vector<char> active(Naxis * nt);

    for (int i = 0; i < N; i++) {
        double ss = 0.0;
        for (int j = 0; j < Naxis; j++) {
            th[j] = stim_thresh_vec(j * N + i);
            double S = 0.0;
            for (int k = 0; k < nt; k++) {
                const double v = X(j * nt * N + i + k * N);
                b[j * nt + k] = std::abs(v);
                sgn[j * nt + k] = (v < 0.0) ? -1.0 : 1.0;
                S += std::abs(v);
            }
            const double u = S / th[j];
            ss += u * u;
        }
        if (sqrt(ss) <= upper) continue;

        double lam = 0.0;
        for (int it = 0; it < 30; it++) {
            // active set per axis at this lam (at most nt drops, since each round removes one or more)
            for (int j = 0; j < Naxis; j++) {
                for (int k = 0; k < nt; k++) active[j * nt + k] = 1;
                double Bj = 0.0;
                int nj = 0;
                for (int k = 0; k < nt; k++) {
                    Bj += b[j * nt + k];
                    nj++;
                }
                for (int round = 0; round < nt; round++) {
                    const double den = th[j] * th[j] + lam * nj;
                    const double delta = (den > 0.0) ? lam * Bj / den : 0.0;
                    bool dropped = false;
                    for (int k = 0; k < nt; k++) {
                        if (active[j * nt + k] && b[j * nt + k] <= delta) {
                            active[j * nt + k] = 0;
                            Bj -= b[j * nt + k];
                            nj--;
                            dropped = true;
                        }
                    }
                    if (!dropped || nj == 0) break;
                }
                B[j] = Bj;
                nact[j] = nj;
            }

            double f = -upper * upper, df = 0.0;
            for (int j = 0; j < Naxis; j++) {
                const double t2 = th[j] * th[j];
                const double den = t2 + lam * nact[j];
                if (den <= 0.0) continue;
                const double s = B[j] * t2 / den;
                f += (s / th[j]) * (s / th[j]);
                df += -2.0 * nact[j] * B[j] * B[j] * t2 / (den * den * den);
            }
            if (df >= 0.0) break; // f is decreasing in lam; guard against round-off
            const double next = lam - f / df;
            const bool done = std::abs(next - lam) <= 1e-15 * (1.0 + std::abs(next));
            lam = (next > 0.0) ? next : 0.5 * lam;
            if (done || std::abs(f) <= 1e-15 * upper * upper) break;
        }

        for (int j = 0; j < Naxis; j++) {
            const double den = th[j] * th[j] + lam * nact[j];
            const double delta = (den > 0.0) ? lam * B[j] / den : 0.0;
            for (int k = 0; k < nt; k++) {
                const double w = std::max(b[j * nt + k] - delta, 0.0);
                X(j * nt * N + i + k * N) = sgn[j * nt + k] * w;
            }
        }
    }

    if (do_equil) {
        X.array() *= eq_rows.array();
    }
    X.array() /= spec_norm;
}

void Op_SAFE::prox(Eigen::VectorXd &X) {
    spdlog::trace("Starting Op_SAFE::prox");

    if (signed_terms13) { // terms 1 and 3 arrive signed; their |.| is taken in the projection
        prox_signed(X);
        return;
    }

    if (do_equil) {
        X.array() /= eq_rows.array();
    }
    X.array() *= spec_norm;

    // METRIC projection onto the combined limit at each sample -- the nearest point, not a rescaling.
    //
    // forward() has already folded a1..a3, g_scale and the |.| into the n_terms blocks, so the constraint
    // sees each axis only through m_j = sum_b x_jb: the set is {x : sum_j (m_j / thresh_j)^2 <= upper^2},
    // the preimage of a ball under a linear map. Only m may change, and since all n_terms blocks enter m
    // with weight 1, the least-norm way to change it is an EQUAL SHARE to each block. Scaling every block
    // by a common factor also lands on the boundary, but it moves the iterate ~20% farther in a direction
    // ~30 deg off the projection, and ADMM's convergence rests on this being the true prox.
    //
    // m_j = m0_j / (1 + lam / thresh_j^2) with lam >= 0 from the KKT conditions; solve
    // sum_j (m_j/thresh_j)^2 = upper^2 by Newton (one step when the thresholds are equal, the usual case).
    const double upper = 1.0 - cushion; // in units of "fraction of the limit"

    std::vector<double> m(Naxis), th(Naxis);
    for (int i = 0; i < N; i++) {
        double ss = 0.0;
        for (int j = 0; j < Naxis; j++) {
            m[j] = axis_stim(X, j, i);
            th[j] = stim_thresh_vec(j * N + i);
            double u = m[j] / th[j];
            ss += u * u;
        }
        if (sqrt(ss) <= upper) continue;

        double lam = 0.0;
        for (int it = 0; it < 50; it++) {
            double f = -upper * upper, df = 0.0;
            for (int j = 0; j < Naxis; j++) {
                const double t2 = th[j] * th[j];
                const double d = 1.0 + lam / t2;
                const double a = m[j] / (th[j] * d);
                f += a * a;
                df += -2.0 * a * a / (d * t2);
            }
            if (df >= 0.0) break; // f is strictly decreasing in lam; guard against round-off
            const double next = lam - f / df;
            const bool done = std::abs(next - lam) <= 1e-15 * (1.0 + std::abs(next));
            lam = (next > 0.0) ? next : 0.5 * lam;
            if (done || std::abs(f) <= 1e-15 * upper * upper) break;
        }
        for (int j = 0; j < Naxis; j++) {
            const double share = (m[j] / (1.0 + lam / (th[j] * th[j])) - m[j]) / n_terms;
            for (int b = 0; b < n_terms; b++) X(j * n_terms * N + i + b * N) += share;
        }
    }

    if (do_equil) {
        X.array() *= eq_rows.array();
    }
    X.array() /= spec_norm;

    spdlog::trace("Finished Op_SAFE::prox");
}

void Op_SAFE::check(Eigen::VectorXd &X) {
    int is_feas = 1;

    if (do_equil) {
        X.array() /= eq_rows.array();
    }
    X.array() *= spec_norm;

    for (int i = 0; i < N; i++) {
        double ss = 0.0;
        for (int j = 0; j < Naxis; j++) {
            double u = axis_stim(X, j, i) / stim_thresh_vec(j * N + i);
            ss += u * u;
        }
        if (sqrt(ss) > 1.0) is_feas = 0;
    }

    if (do_equil) {
        X.array() *= eq_rows.array();
    }
    X.array() /= spec_norm;

    hist_feas.push_back(is_feas);
}

void Op_SAFE::freeze_linearization(Eigen::VectorXd &X) {
    // Capture the true |.| signs at X, then hold them so the CG sees a fixed linear (symmetric) system;
    // unfreeze_linearization() restores the true forward for prox and check.
    freeze_signs = false;
    Ax_temp.setZero(Ax_size);
    forward_op(X, Ax_temp);
    freeze_signs = true;
}

double Op_SAFE::linearization_error(const Eigen::VectorXd &x_new) {
    // Relative mismatch between the frozen-linear forward at x_new (what the CG assumed) and the true
    // nonlinear forward. Leaves the frozen signs and freeze_signs unchanged.
    Eigen::VectorXd xc = x_new;
    Eigen::VectorXd pred(Ax_size), act(Ax_size);
    const bool saved = freeze_signs;
    const Eigen::VectorXd s1 = signs1, s2 = signs2, s3 = signs3;
    freeze_signs = true;
    forward_op(xc, pred); // frozen-linear prediction
    freeze_signs = false;
    forward_op(xc, act);  // true nonlinear SAFE
    signs1 = s1;
    signs2 = s2;
    signs3 = s3;
    freeze_signs = saved;
    double an = act.norm();
    return (an > 0.0) ? (act - pred).norm() / an : 0.0;
}

double Op_SAFE::constraint_violation(const Eigen::VectorXd &x_new) {
    // max over samples of the fractional excess over the combined limit, sqrt(sum_j (stim_j/thresh_j)^2) - 1,
    // matching check(); leaves the frozen-sign state unchanged.
    Eigen::VectorXd xc = x_new;
    Eigen::VectorXd out(Ax_size);
    const bool saved = freeze_signs;
    const Eigen::VectorXd s1 = signs1, s2 = signs2, s3 = signs3;
    freeze_signs = false;
    forward(xc, out);
    signs1 = s1;
    signs2 = s2;
    signs3 = s3;
    freeze_signs = saved;
    double viol = 0.0;
    for (int i = 0; i < N; i++) {
        double ss = 0.0;
        for (int j = 0; j < Naxis; j++) {
            double u = axis_stim(out, j, i) / stim_thresh_vec(j * N + i);
            ss += u * u;
        }
        double over = sqrt(ss) - 1.0;
        if (over > viol) viol = over;
    }
    return viol;
}

void SAFEParams::set_demo_params() {
    tau1[0] = 0.14 / 1000.0;
    tau2[0] = 12.0 / 1000.0;
    tau3[0] = 0.52 / 1000.0;
    a1[0] = 0.32;
    a2[0] = 0.20;
    a3[0] = 0.48;
    stim_limit[0] = 28;
    g_scale[0] = 0.34;

    tau1[1] = 1.5 / 1000.0;
    tau2[1] = 2.5 / 1000.0;
    tau3[1] = 0.15 / 1000.0;
    a1[1] = 0.55;
    a2[1] = 0.15;
    a3[1] = 0.3;
    stim_limit[1] = 15;
    g_scale[1] = 0.31;

    tau1[2] = 2 / 1000.0;
    tau2[2] = 0.12 / 1000.0;
    tau3[2] = 1 / 1000.0;
    a1[2] = 0.42;
    a2[2] = 0.4;
    a3[2] = 0.18;
    stim_limit[2] = 25;
    g_scale[2] = 0.25;
}

void SAFEParams::set_params(const Eigen::VectorXd &_tau1, const Eigen::VectorXd &_tau2, const Eigen::VectorXd &_tau3,
                            const Eigen::VectorXd &_a1, const Eigen::VectorXd &_a2, const Eigen::VectorXd &_a3,
                            const Eigen::VectorXd &_stim_limit, const Eigen::VectorXd &_g_scale) {
    for (int i = 0; i < 3; i++) {
        tau1[i] = _tau1(i);
        tau2[i] = _tau2(i);
        tau3[i] = _tau3(i);
        a1[i] = _a1(i);
        a2[i] = _a2(i);
        a3[i] = _a3(i);
        stim_limit[i] = _stim_limit(i);
        g_scale[i] = _g_scale(i);
    }
}

void SAFEParams::swap_first_axes(int new_first_axis) {
    if (new_first_axis < 0 || new_first_axis >= 3) {
        spdlog::error("Invalid axis index: {}", new_first_axis);
        return;
    } else if (new_first_axis > 0) {
        std::swap(tau1[0], tau1[new_first_axis]);
        std::swap(tau2[0], tau2[new_first_axis]);
        std::swap(tau3[0], tau3[new_first_axis]);
        std::swap(a1[0], a1[new_first_axis]);
        std::swap(a2[0], a2[new_first_axis]);
        std::swap(a3[0], a3[new_first_axis]);
        std::swap(stim_limit[0], stim_limit[new_first_axis]);
        std::swap(g_scale[0], g_scale[new_first_axis]);
    }
}

void SAFEParams::calc_alphas(double dt) {
    for (int i = 0; i < 3; i++) {
        alpha1[i] = dt / (tau1[i] + dt);
        alpha2[i] = dt / (tau2[i] + dt);
        alpha3[i] = dt / (tau3[i] + dt);
    }
}

} // namespace Gropt
