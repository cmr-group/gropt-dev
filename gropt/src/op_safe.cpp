#include "spdlog/spdlog.h"

#include "op_safe.hpp"

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

    safe_params.calc_alphas(pdata->dt);

    if ((safe_params.a3[0] == 0) && (safe_params.a3[1] == 0) && (safe_params.a3[2] == 0)) {
        n_terms = 2;
        spdlog::trace("Op_SAFE::init  n_terms = {}", n_terms);
    }

    // We will ignore this for now and just use stim_thresh_vec directly
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

    if (do_init_weights) {
        // weight is set on workspace, not here
    }

    signs1.setZero(pdata->Naxis * pdata->N);
    signs2.setZero(pdata->Naxis * pdata->N);
    signs3.setZero(pdata->Naxis * pdata->N);
    stim1.setZero(pdata->Naxis * pdata->N);
    stim2.setZero(pdata->Naxis * pdata->N);
    stim3.setZero(pdata->Naxis * pdata->N);

    Operator::init();

    // Analytical spec norm was non-trivial, just do the estimate instead.
    // TODO: Test the amount here, or work on stopping early in the spec norm code.
    spec_norm = estimate_self_spec_norm(30);
    spec_norm2 = spec_norm * spec_norm;

    spdlog::trace("Op_SAFE::init  spec_norm = {:.4f}", spec_norm);
    spdlog::trace("Op_SAFE::init  Done!");
}

void Op_SAFE::forward(Eigen::VectorXd &X, Eigen::VectorXd &out) {
    // out = diff(X)/dt
    for (int j = 0; j < Naxis; j++) {
        out(j * N) = X(j * N) / dt;
        for (int i = 1; i < N; i++) {
            out(j * N + i) = (X(j * N + i) - X(j * N + i - 1)) / dt;
        }
    }

    // stim1 = tau_filter_1(dX/dt)
    stim1.setZero();
    for (int j = 0; j < Naxis; j++) {
        stim1(j * N) = safe_params.alpha1[j] * out(j * N);
        for (int i = 1; i < N; i++) {
            stim1(j * N + i) =
                safe_params.alpha1[j] * out(j * N + i) + (1.0 - safe_params.alpha1[j]) * stim1(j * N + i - 1);
        }
    }

    // stim1 = abs(tau_filter_1(dX/dt))   (softabs + smooth sign when safe_eps>0)
    for (int i = 0; i < stim1.size(); i++) {
        double v = stim1(i);
        if (freeze_signs) {
            stim1(i) = signs1(i) * v; // frozen linearization: apply held sign LINEARLY (no recapture)
        } else if (safe_eps > 0.0) {
            double sa = sqrt(v * v + safe_eps * safe_eps);
            signs1(i) = v / sa; // smooth sign in [-1,1], continuous through 0 (= d/dv of softabs)
            stim1(i) = sa;
        } else {
            signs1(i) = (v < 0.0) ? -1.0 : 1.0;
            stim1(i) = (v < 0.0) ? -v : v;
        }
    }

    // stim2 = dX/dt
    for (int j = 0; j < Naxis; j++) {
        for (int i = 0; i < N; i++) {
            stim2(j * N + i) = out(j * N + i);
        }
    }

    // stim2 = abs(dX/dt)   (softabs + smooth sign when safe_eps>0)
    for (int i = 0; i < stim1.size(); i++) {
        double v = stim2(i);
        if (freeze_signs) {
            stim2(i) = signs2(i) * v; // frozen linearization (path 2: sign is applied BEFORE filter2)
        } else if (safe_eps > 0.0) {
            double sa = sqrt(v * v + safe_eps * safe_eps);
            signs2(i) = v / sa;
            stim2(i) = sa;
        } else {
            signs2(i) = (v < 0.0) ? -1.0 : 1.0;
            stim2(i) = (v < 0.0) ? -v : v;
        }
    }

    // stim2 = tau_filter_2(abs(dX/dt))
    for (int j = 0; j < Naxis; j++) {
        stim2(j * N) = safe_params.alpha2[j] * stim2(j * N);
        for (int i = 1; i < N; i++) {
            stim2(j * N + i) =
                safe_params.alpha2[j] * stim2(j * N + i) + (1.0 - safe_params.alpha2[j]) * stim2(j * N + i - 1);
        }
    }

    // stim3 = tau_filter_3(dX/dt)
    stim3.setZero();
    for (int j = 0; j < Naxis; j++) {
        stim3(j * N) = safe_params.alpha3[j] * out(j * N);
        for (int i = 1; i < N; i++) {
            stim3(j * N + i) =
                safe_params.alpha3[j] * out(j * N + i) + (1.0 - safe_params.alpha3[j]) * stim3(j * N + i - 1);
        }
    }

    // stim3 = abs(tau_filter_3(dX/dt))   (softabs + smooth sign when safe_eps>0)
    for (int i = 0; i < stim3.size(); i++) {
        double v = stim3(i);
        if (freeze_signs) {
            stim3(i) = signs3(i) * v; // frozen linearization: apply held sign LINEARLY (no recapture)
        } else if (safe_eps > 0.0) {
            double sa = sqrt(v * v + safe_eps * safe_eps);
            signs3(i) = v / sa;
            stim3(i) = sa;
        } else {
            signs3(i) = (v < 0.0) ? -1.0 : 1.0;
            stim3(i) = (v < 0.0) ? -v : v;
        }
    }

    // The return vector is [stim1; stim2; stim3] scaled by a, stim_limit, g_scale
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

    // stim1 = signs1 * stim1
    for (int i = 0; i < stim1.size(); i++) {
        stim1(i) = signs1(i) * stim1(i);
    }

    // stim1 = tau_filter_1_T(stim1)
    for (int j = 0; j < Naxis; j++) {
        stim1(j * N + N - 1) = safe_params.alpha1[j] * stim1(j * N + N - 1);
        for (int i = N - 2; i >= 0; i--) {
            stim1(j * N + i) =
                safe_params.alpha1[j] * stim1(j * N + i) + (1 - safe_params.alpha1[j]) * stim1(j * N + i + 1);
        }
    }

    // stim2 = tau_filter_2_T(stim2)
    for (int j = 0; j < Naxis; j++) {
        stim2(j * N + N - 1) = safe_params.alpha2[j] * stim2(j * N + N - 1);
        for (int i = N - 2; i >= 0; i--) {
            stim2(j * N + i) =
                safe_params.alpha2[j] * stim2(j * N + i) + (1 - safe_params.alpha2[j]) * stim2(j * N + i + 1);
        }
    }

    // stim2 = signs2 * stim2
    for (int i = 0; i < stim2.size(); i++) {
        stim2(i) = signs2(i) * stim2(i);
    }

    if (n_terms == 3) {

        // stim3 = signs3 * stim3
        for (int i = 0; i < stim3.size(); i++) {
            stim3(i) = signs3(i) * stim3(i);
        }

        // stim3 = tau_filter_3_T(stim3)
        for (int j = 0; j < Naxis; j++) {
            stim3(j * N + N - 1) = safe_params.alpha3[j] * stim3(j * N + N - 1);
            for (int i = N - 2; i >= 0; i--) {
                stim3(j * N + i) =
                    safe_params.alpha3[j] * stim3(j * N + i) + (1 - safe_params.alpha3[j]) * stim3(j * N + i + 1);
            }
        }
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

    // out = diff(stim1 + stim2 + stim3)/dt
    for (int j = 0; j < Naxis; j++) {
        for (int i = 0; i < N - 1; i++) {
            out(j * N + i) = (out(j * N + i) - out(j * N + i + 1)) / dt;
        }
        out(j * N + N - 1) = out(j * N + N - 1) / dt;
    }
}

void Op_SAFE::prox(Eigen::VectorXd &X) {
    spdlog::trace("Starting Op_SAFE::prox");

    if (do_equil) {
        X.array() /= eq_rows.array();
    }
    X.array() *= spec_norm;

    // This just sums the three terms [stim1, stim2, stim3], it is not the cross axis sum.
    x_temp.setZero();
    for (int j = 0; j < Naxis; j++) {
        for (int i = 0; i < N; i++) {
            if (n_terms == 3) {
                x_temp(j * N + i) =
                    X(j * n_terms * N + i) + X(j * n_terms * N + i + N) + X(j * n_terms * N + i + 2 * N);
            } else if (n_terms == 2) {
                x_temp(j * N + i) = X(j * n_terms * N + i) + X(j * n_terms * N + i + N);
            }
        }
    }

    if (rot_variant) {
        for (int j = 0; j < Naxis; j++) {
            for (int i = 0; i < N; i++) {
                double val = abs(x_temp(j * N + i));

                double upper_bound = (1 - cushion) * stim_thresh_vec(j * N + i);
                if (val > upper_bound) {
                    X(j * n_terms * N + i) *= (upper_bound / val);
                    X(j * n_terms * N + i + N) *= (upper_bound / val);
                    if (n_terms == 3) {
                        X(j * n_terms * N + i + 2 * N) *= (upper_bound / val);
                    }
                }
            }
        }
    } else {
        for (int i = 0; i < N; i++) {
            double val = 0.0;
            for (int j = 0; j < Naxis; j++) {
                val += X(j * N + i) * X(j * N + i);
            }
            val = sqrt(val);
            double upper_bound = (1 - cushion) * stim_thresh_vec(i);

            if (val > upper_bound) {
                for (int j = 0; j < Naxis; j++) {
                    X(j * N + i) *= (upper_bound / val);
                }
            }
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

    x_temp.setZero();
    for (int j = 0; j < Naxis; j++) {
        for (int i = 0; i < N; i++) {
            if (n_terms == 3) {
                x_temp(j * N + i) =
                    X(j * n_terms * N + i) + X(j * n_terms * N + i + N) + X(j * n_terms * N + i + 2 * N);
            } else if (n_terms == 2) {
                x_temp(j * N + i) = X(j * n_terms * N + i) + X(j * n_terms * N + i + N);
            }
        }
    }

    if (rot_variant) {
        for (int j = 0; j < Naxis; j++) {
            for (int i = 0; i < N; i++) {
                double val = abs(x_temp(j * N + i));
                double upper_bound = stim_thresh_vec(j * N + i);

                if (val > upper_bound) {
                    is_feas = 0;
                }
            }
        }
    } else {
        for (int i = 0; i < N; i++) {
            double val = 0.0;
            for (int j = 0; j < Naxis; j++) {
                val += X(j * N + i) * X(j * N + i);
            }
            val = sqrt(val);
            double upper_bound = stim_thresh_vec(i);

            if (val > upper_bound) {
                is_feas = 0;
            }
        }
    }

    if (do_equil) {
        X.array() *= eq_rows.array();
    }
    X.array() /= spec_norm;

    hist_feas.push_back(is_feas);
}

void Op_SAFE::freeze_linearization(Eigen::VectorXd &X) {
    // Recapture the true |.| signs once at the current outer iterate X (nonlinear forward), then hold
    // them fixed. During the inner CG, forward() then applies those signs LINEARLY, so the CG's LHS is a
    // consistent symmetric operator instead of re-capturing signs off every search direction. Restored by
    // unfreeze_linearization() after the solve, so the prox / feasibility check use the true forward.
    freeze_signs = false;
    Ax_temp.setZero(Ax_size);
    forward_op(X, Ax_temp);
    freeze_signs = true;
}

double Op_SAFE::linearization_error(const Eigen::VectorXd &x_new) {
    // signs1/2/3 currently hold the linearization frozen at the outer iterate. Compare the FROZEN-LINEAR
    // forward at x_new (what the CG's clamp assumed) against the TRUE nonlinear forward at x_new. Relative
    // mismatch = how far the step walked out of the linearization's valid region. NOTE: the true forward
    // recaptures signs (overwriting the frozen ones), which is fine -- the caller unfreezes / re-freezes
    // right after. Evaluate the LINEAR prediction FIRST while the frozen signs are still intact.
    Eigen::VectorXd xc = x_new;
    Eigen::VectorXd pred(Ax_size), act(Ax_size);
    bool saved = freeze_signs;
    freeze_signs = true;
    forward_op(xc, pred); // frozen-linear prediction (signs held at the outer iterate)
    freeze_signs = false;
    forward_op(xc, act); // true nonlinear SAFE (recaptures signs)
    freeze_signs = saved;
    double an = act.norm();
    return (an > 0.0) ? (act - pred).norm() / an : 0.0;
}

double Op_SAFE::constraint_violation(const Eigen::VectorXd &x_new) {
    // WORST-sample TRUE (nonlinear) SAFE overage: max over samples of (|SAFE(x)| - limit), 0 if all
    // feasible. Scale-free (SAFE units, ~[0,limit]) so a single feasibility slack works across problems.
    // Uses the raw forward (physical scaled terms, NOT /spec_norm) with the true |.|, matching check().
    Eigen::VectorXd xc = x_new;
    Eigen::VectorXd out(Ax_size);
    bool saved = freeze_signs;
    freeze_signs = false;
    forward(xc, out);
    freeze_signs = saved;
    double viol = 0.0;
    for (int j = 0; j < Naxis; j++) {
        for (int i = 0; i < N; i++) {
            double val;
            if (n_terms == 3) {
                val = out(j * n_terms * N + i) + out(j * n_terms * N + i + N) + out(j * n_terms * N + i + 2 * N);
            } else {
                val = out(j * n_terms * N + i) + out(j * n_terms * N + i + N);
            }
            double av = (val < 0.0) ? -val : val;
            double over = av - stim_thresh_vec(j * N + i);
            if (over > viol) viol = over;
        }
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
