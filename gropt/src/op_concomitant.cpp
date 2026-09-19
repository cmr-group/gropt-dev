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

    Operator::init();
}

void Op_Concomitant::append_eq_rows(std::vector<Eigen::VectorXd> &rows, std::vector<double> &targets,
                                    const Eigen::VectorXd &x0) const {
    if (!use_projection) return;

    // c is quadratic, so a·x0 = 2 c(x0) and the linearization c(x0) + a·(x - x0) = 0 becomes a·x = c(x0).
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

    // Skip while the waveform is still essentially zero
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

} // namespace Gropt
