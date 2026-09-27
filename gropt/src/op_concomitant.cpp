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

// pos/neg use disjoint samples so prox can rescale one side alone; exact_quad drops the interval
// straddling the flip, which lies inside the zeroed 180 window.
void Op_Concomitant::energies(const Eigen::VectorXd &X, double &pos, double &neg, Eigen::VectorXd *d_pos,
                              Eigen::VectorXd *d_neg) const {
    pos = 0.0;
    neg = 0.0;
    if (d_pos) d_pos->setZero(X.size());
    if (d_neg) d_neg->setZero(X.size());

    if (!exact_quad) { // legacy: rectangle sum over the flattened vector
        for (int i = start_idx; i < X.size(); i++) {
            if (pdata->inv_vec(i) > 0) {
                pos += X(i) * X(i) * pdata->dt;
                if (d_pos) (*d_pos)(i) = 2.0 * pdata->dt * X(i);
            } else if (pdata->inv_vec(i) < 0) {
                neg += X(i) * X(i) * pdata->dt;
                if (d_neg) (*d_neg)(i) = 2.0 * pdata->dt * X(i);
            }
        }
        return;
    }

    const int N = pdata->N;
    const double dt = pdata->dt;
    for (int ax = 0; ax < pdata->Naxis; ax++) {
        const int off = ax * N;
        const int i0 = (start_idx > off) ? start_idx - off : 0;
        for (int i = i0; i < N - 1; i++) {
            const double a = X(off + i);
            const double b = X(off + i + 1);
            const bool is_pos = pdata->inv_vec(off + i) > 0 && pdata->inv_vec(off + i + 1) > 0;
            if (!is_pos && !(pdata->inv_vec(off + i) < 0 && pdata->inv_vec(off + i + 1) < 0)) continue;
            const double seg = dt / 3.0 * (a * a + a * b + b * b);
            Eigen::VectorXd *d = is_pos ? d_pos : d_neg;
            (is_pos ? pos : neg) += seg;
            if (d) {
                (*d)(off + i) += dt / 3.0 * (2.0 * a + b);
                (*d)(off + i + 1) += dt / 3.0 * (a + 2.0 * b);
            }
        }
    }
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
    double pos = 0.0, neg = 0.0;
    Eigen::VectorXd d_pos, d_neg;
    energies(x0, pos, neg, &d_pos, &d_neg);
    Eigen::VectorXd a = d_pos - target * d_neg;

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

    double pos = 0.0, neg = 0.0;
    energies(X, pos, neg, nullptr, nullptr);

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

    double pos = 0.0, neg = 0.0;
    energies(X, pos, neg, nullptr, nullptr);

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
