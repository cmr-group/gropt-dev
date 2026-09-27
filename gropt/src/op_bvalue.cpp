#include "spdlog/spdlog.h"

#include "op_bvalue.hpp"

namespace Gropt {

Op_BValue::Op_BValue(const ProblemData &_pdata, double _bval_target, double _bval_tol0, int _start_idx0, int _stop_idx0,
                     double _weight_mod, BVALUE_MODE _mode, double _max_scale)
    : Operator(_pdata) {
    name = "b-value";

    bval_target = _bval_target;
    bval_tol0 = _bval_tol0;

    start_idx0 = _start_idx0;
    stop_idx0 = _stop_idx0;

    start_idx = start_idx0;
    stop_idx = stop_idx0;

    weight_mod = _weight_mod;

    mode = _mode;
    max_scale = _max_scale;
}

void Op_BValue::init() {
    spdlog::trace("Op_BValue::init  N = {}", pdata->N);

    target = bval_target;
    tol0 = bval_tol0;
    tol = (1.0 - cushion) * tol0;

    GAMMA = 267.5221900e6;                                          // [rad/s/T]
    MAT_SCALE = pow((GAMMA / 1000.0 * pdata->dt), 2.0) * pdata->dt; // 1/1000 is for m->mm in b-value

    // start/stop <= 0: whole axis
    if (start_idx <= 0) {
        i_start = 0;
    } else {
        i_start = start_idx;
    }

    if (stop_idx <= 0) {
        i_stop = pdata->N;
    } else {
        i_stop = stop_idx;
    }

    // Legendre scales, so that b = ||Ax||^2 with out = [A0 c0, A1 c1, A2 c2] per interval
    const double gam = GAMMA / 1000.0;
    PWL_A0 = gam * sqrt(pdata->dt);
    PWL_A1 = gam * sqrt(pdata->dt / 3.0);
    PWL_A2 = gam * sqrt(pdata->dt / 5.0);

    ax_per_axis = pwl_quad ? 3 * (pdata->N - 1) : pdata->N;
    Ax_size = pdata->Naxis * ax_per_axis;

    if (!pwl_quad) {
        int Nnorm = i_stop - i_start;
        spec_norm2 = (Nnorm * Nnorm + Nnorm) / 2.0 * MAT_SCALE * 0.1175 * 4;
        spec_norm = sqrt(spec_norm2);
    }

    if (do_init_weights) {
        obj_weight = -1.0;
        obj_weight *= weight_mod;
    }

    Operator::init();

    if (pwl_quad) {
        // the closed form above is for the cumsum operator; this A is a different shape
        spec_norm = estimate_self_spec_norm(30);
        spec_norm2 = spec_norm * spec_norm;
        spdlog::trace("Op_BValue::init  pwl_quad spec_norm = {:.4e}", spec_norm);
    }
}

void Op_BValue::forward(Eigen::VectorXd &X, Eigen::VectorXd &out) {
    out.setZero();

    if (pwl_quad) {
        const double h = pdata->dt;
        for (int j = 0; j < Naxis; j++) {
            const int jN = j * N, base = j * ax_per_axis;
            double Q = 0.0; // integral of g up to node i, by the trapezoid rule (exact here)
            for (int i = i_start; i <= i_stop - 2; i++) {
                const double gi = X(jN + i) * pdata->inv_vec(jN + i);
                const double gn = X(jN + i + 1) * pdata->inv_vec(jN + i + 1);
                out(base + 3 * i + 0) = PWL_A0 * (Q + h * (2.0 * gi + gn) / 6.0);
                out(base + 3 * i + 1) = PWL_A1 * (h * (gi + gn) / 4.0);
                out(base + 3 * i + 2) = PWL_A2 * (h * (gn - gi) / 12.0);
                Q += h * (gi + gn) / 2.0;
            }
        }
        return;
    }

    for (int j = 0; j < Naxis; j++) {
        int jN = j * N;
        double gt = 0;
        for (int i = i_start; i < i_stop; i++) {
            gt += X(jN + i) * pdata->inv_vec(jN + i);
            out(jN + i) = gt * sqrt(MAT_SCALE);
        }
    }
}

void Op_BValue::transpose(Eigen::VectorXd &X, Eigen::VectorXd &out) {
    out.setZero();

    if (pwl_quad) {
        // Adjoint of the forward above. Q_i is a trapezoid partial sum, so its transpose is a suffix
        // sum of the c0 duals (half weight at the window ends); c1/c2 are local to the two nodes of
        // their own interval.
        const double h = pdata->dt;
        for (int j = 0; j < Naxis; j++) {
            const int jN = j * N, base = j * ax_per_axis;
            const int i_last = i_stop - 2; // last interval
            if (i_last < i_start) continue;
            double suf = 0.0; // sum of c0 duals over intervals strictly after the current node
            for (int k = i_stop - 1; k >= i_start; k--) {
                const double u_k = (k <= i_last) ? PWL_A0 * X(base + 3 * k) : 0.0;
                double acc = (k == i_start) ? h * 0.5 * suf : h * (0.5 * u_k + suf);

                if (k <= i_last) { // node k is the left end of interval k
                    acc += u_k * h / 3.0 + PWL_A1 * X(base + 3 * k + 1) * h / 4.0 -
                           PWL_A2 * X(base + 3 * k + 2) * h / 12.0;
                }
                if (k - 1 >= i_start) { // node k is the right end of interval k-1
                    acc += PWL_A0 * X(base + 3 * (k - 1)) * h / 6.0 +
                           PWL_A1 * X(base + 3 * (k - 1) + 1) * h / 4.0 +
                           PWL_A2 * X(base + 3 * (k - 1) + 2) * h / 12.0;
                }
                out(jN + k) = acc * pdata->inv_vec(jN + k);
                if (k <= i_last) suf += u_k;
            }
        }
        return;
    }

    for (int j = 0; j < Naxis; j++) {
        int jN = j * N;
        double gt = 0;
        for (int i = i_stop - 1; i >= i_start; i--) {
            gt += X(jN + i) * sqrt(MAT_SCALE);
            out(jN + i) = gt * pdata->inv_vec(jN + i);
        }
    }
}

// TODO: This seems all wrong for three axis case, need to fix
void Op_BValue::prox(Eigen::VectorXd &X) {
    spdlog::trace("Starting Op_BValue::prox");

    if (do_equil) {
        X.array() /= eq_rows.array();
    }
    X.array() *= spec_norm;

    for (int j = 0; j < Naxis; j++) {
        // b for this axis is the squared norm of its Ax block, whichever quadrature built it, so the
        // radial rescaling below is the same projection in both cases
        double xnorm = X.segment(j * ax_per_axis, ax_per_axis).norm();

        if (mode == BVALUE_MODE_MINVAL) {
            double min_val = sqrt(target);

            if (xnorm < min_val) {
                X.segment(j * ax_per_axis, ax_per_axis) *= (min_val / xnorm);
            }
        } else if (mode == BVALUE_MODE_MINVALMAX) {
            double min_val = sqrt(target);

            if (xnorm < min_val) {
                X.segment(j * ax_per_axis, ax_per_axis) *= (max_scale * min_val / xnorm);
            } else {
                X.segment(j * ax_per_axis, ax_per_axis) *= max_scale;
            }
        } else if (mode == BVALUE_MODE_SETVAL) {
            double min_val = sqrt(target - tol);
            double max_val = sqrt(target + tol);

            if (xnorm < min_val) {
                X.segment(j * ax_per_axis, ax_per_axis) *= (min_val / xnorm);
            } else if (xnorm > max_val) {
                X.segment(j * ax_per_axis, ax_per_axis) *= (max_val / xnorm);
            }
        } else {
            spdlog::error("Unknown BVALUE_MODE in Op_BValue::prox");
        }
    }

    if (do_equil) {
        X.array() *= eq_rows.array();
    }
    X.array() /= spec_norm;

    spdlog::trace("Finished Op_BValue::prox");
}

void Op_BValue::check(Eigen::VectorXd &X) {
    int is_feas = 1;

    if (do_equil) {
        X.array() /= eq_rows.array();
    }
    X.array() *= spec_norm;

    for (int j = 0; j < Naxis; j++) {
        double bval_t = (X.segment(j * ax_per_axis, ax_per_axis)).squaredNorm();

        if (mode == BVALUE_MODE_MINVAL || mode == BVALUE_MODE_MINVALMAX) {
            if (bval_t < target) {
                is_feas = 0;
            }
        } else if (mode == BVALUE_MODE_SETVAL) {
            double d_bval = fabs(bval_t - target);
            if (d_bval > tol0) {
                is_feas = 0;
            }
        } else {
            spdlog::error("Unknown BVALUE_MODE in Op_BValue::check");
        }
    }

    if (do_equil) {
        X.array() *= eq_rows.array();
    }
    X.array() /= spec_norm;

    hist_feas.push_back(is_feas);
}

double Op_BValue::get_bvalue(Eigen::VectorXd &X) {
    Ax_temp.setZero();
    forward_op(X, Ax_temp);
    Ax_temp.array() *= spec_norm;

    return Ax_temp.squaredNorm();
}

} // namespace Gropt
