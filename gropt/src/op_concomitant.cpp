#include <cmath>

#include "spdlog/spdlog.h"

#include "op_concomitant.hpp"

namespace Gropt {

Op_Concomitant::Op_Concomitant(const ProblemData &_pdata, int _start_idx, bool _rot_variant, double _weight_mod,
                               double _tol0, double _target)
    : Operator(_pdata) {
    name = "Concomitant";

    start_idx = _start_idx;
    rot_variant = _rot_variant;
    weight_mod = _weight_mod;
    con_tol0 = _tol0;
    con_target = _target;
}

void Op_Concomitant::init() {
    spdlog::trace("Op_Concomitant::init  N = {}", pdata->N);

    target = con_target;
    tol0 = con_tol0;
    tol = (1.0 - cushion) * tol0;

    spec_norm2 = 1.0;
    spec_norm = 1.0;

    Ax_size = pdata->Naxis * pdata->N;

    if (do_init_weights) {
        obj_weight = 1.0;
        obj_weight *= weight_mod;
    }

    // Augmented LAgrangian worked, but never really helped...
    auglag_rho = weight_mod; // objective-mode penalty parameter (only used when added to all_obj)
    auglag_lambda = 0.0;

    Operator::init();
}

void Op_Concomitant::append_eq_rows(std::vector<Eigen::VectorXd> &rows, std::vector<double> &targets,
                                    const Eigen::VectorXd &x0) const {
    if (!use_projection) return;

    // Linearize c(x) = pos - target·neg = (Σ_pre x²dt − target·Σ_post x²dt) at x0 (SQP step); c=0 is
    // pos/neg = target. The linear row is a = grad c(x0) = 2·dt·x0 on pre, −2·target·dt·x0 on post;
    // the row target works out to c0 = pos0 − target·neg0 (since a·x0 = 2·c0). As x0 -> feasible,
    // c0 -> 0 and the row enforces the tangent grad c·x = c0.
    Eigen::VectorXd a = Eigen::VectorXd::Zero(x0.size());
    double pos = 0.0;
    double neg = 0.0;
    for (int i = start_idx; i < x0.size(); i++) {
        if (pdata->inv_vec(i) > 0) {
            pos += x0(i) * x0(i) * pdata->dt;
            a(i) = 2.0 * pdata->dt * x0(i);
        } else if (pdata->inv_vec(i) < 0) {
            neg += x0(i) * x0(i) * pdata->dt;
            a(i) = -2.0 * target * pdata->dt * x0(i);
        }
    }

    // Skip while the waveform is still essentially zero.
    if (a.norm() < 1e-12) return;

    rows.push_back(a);
    targets.push_back(pos - target * neg);
}

void Op_Concomitant::forward(Eigen::VectorXd &X, Eigen::VectorXd &out) { out = X; }

void Op_Concomitant::transpose(Eigen::VectorXd &X, Eigen::VectorXd &out) { out = X; }

void Op_Concomitant::prox(Eigen::VectorXd &X) {
    spdlog::trace("Starting Op_Concomitant::prox");

    if (do_equil) {
        X.array() /= eq_rows.array();
    }
    X.array() *= spec_norm;

    double pos = 0.0;
    double neg = 0.0;

    for (int i = start_idx; i < X.size(); i++) {
        if (pdata->inv_vec(i) > 0) {
            pos += X(i) * X(i) * pdata->dt;
        } else if (pdata->inv_vec(i) < 0) {
            neg += X(i) * X(i) * pdata->dt;
        }
    }

    double eps = 1e-12 * (pos + neg);
    if (pos > eps && neg > eps) {
        double q = pos / neg;
        double lo = target - tol;
        double hi = target + tol;
        double q_clamp = (q < lo) ? lo : ((q > hi) ? hi : q);
        if (q_clamp != q) {
            double a = std::sqrt(pos); // ||u|| (pre)
            double b = std::sqrt(neg); // ||v|| (post)
            double c = std::sqrt(q_clamp); // target norm ratio ||u'||/||v'||
            double rv = (c * a + b) / (c * c + 1.0); // new post norm (minimal displacement)
            double ru = c * rv;                      // new pre norm
            double s_pos = ru / a;
            double s_neg = rv / b;
            for (int i = start_idx; i < X.size(); i++) {
                if (pdata->inv_vec(i) > 0) {
                    X(i) *= s_pos;
                } else if (pdata->inv_vec(i) < 0) {
                    X(i) *= s_neg;
                }
            }
        }
    }

    if (do_equil) {
        X.array() *= eq_rows.array();
    }
    X.array() /= spec_norm;

    spdlog::trace("Finished Op_Concomitant::prox");
}

void Op_Concomitant::check(Eigen::VectorXd &X) {
    int is_feas = 1;

    if (do_equil) {
        X.array() /= eq_rows.array();
    }
    X.array() *= spec_norm;

    double pos = 0.0;
    double neg = 0.0;

    for (int i = start_idx; i < X.size(); i++) {
        if (pdata->inv_vec(i) > 0) {
            pos += X(i) * X(i) * pdata->dt;
        } else if (pdata->inv_vec(i) < 0) {
            neg += X(i) * X(i) * pdata->dt;
        }
    }

    double eps = 1e-12 * (pos + neg);
    if (pos <= eps || neg <= eps) {
        is_feas = 0; // one side ~0: degenerate, not a balanced state
    } else {
        double q = pos / neg;
        if ((q < target - tol0) || (q > target + tol0)) {
            is_feas = 0;
        }
    }

    if (do_equil) {
        X.array() *= eq_rows.array();
    }
    X.array() /= spec_norm;

    hist_feas.push_back(is_feas);
}

// ---- Objective (all_obj) mode: augmented Lagrangian for c(x)=pos-neg=0 ----
// Freeze the linearization at the outer iterate x0 and advance the scalar dual (method of
// multipliers). g is the masked gradient grad c(x0) = 2*dt*x0 (signed), so fixed DOFs aren't pushed.
void Op_Concomitant::update_obj_state(const Eigen::VectorXd &x0) {
    if (auglag_g.size() != x0.size()) auglag_g.setZero(x0.size());
    auglag_g.setZero();

    double pos = 0.0;
    double neg = 0.0;
    for (int i = start_idx; i < x0.size(); i++) {
        if (pdata->inv_vec(i) > 0) {
            pos += x0(i) * x0(i) * pdata->dt;
            auglag_g(i) = 2.0 * pdata->dt * x0(i) * pdata->fixer(i);
        } else if (pdata->inv_vec(i) < 0) {
            neg += x0(i) * x0(i) * pdata->dt;
            auglag_g(i) = -2.0 * target * pdata->dt * x0(i) * pdata->fixer(i);
        }
    }

    auglag_c = pos - target * neg;   // c(x0) = pos - target*neg
    auglag_gx0 = auglag_g.dot(x0);   // g . x0
    auglag_lambda += auglag_rho * auglag_c; // dual update (equality)
}

// Gauss-Newton curvature of (rho/2)c^2 about the frozen linearization: rho * g (g.x). Rank-1 PSD,
// so the CG matrix stays positive-definite.
void Op_Concomitant::add_obj(Eigen::VectorXd &X, Eigen::VectorXd &out) {
    if (auglag_g.size() != X.size()) return;
    out.array() += (auglag_rho * auglag_g.dot(X)) * auglag_g.array();
}

// Frozen pull from the linearized penalty + scalar dual:  (rho*(g.x0 - c0) - lambda) * g.
void Op_Concomitant::add_obj_rhs(Eigen::VectorXd &x0, Eigen::VectorXd &out, bool normalize) {
    if (auglag_g.size() != out.size()) return;
    double coef = auglag_rho * (auglag_gx0 - auglag_c) - auglag_lambda;
    out.array() += coef * auglag_g.array();
}

} // namespace Gropt
