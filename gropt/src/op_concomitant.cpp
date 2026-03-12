#include "spdlog/spdlog.h"

#include "op_concomitant.hpp"

namespace Gropt {

Op_Concomitant::Op_Concomitant(const ProblemData &_pdata, bool _rot_variant, double _weight_mod) : Operator(_pdata) {
    name = "Concomitant";
    rot_variant = _rot_variant;
    weight_mod = _weight_mod;
}

void Op_Concomitant::init() {
    spdlog::trace("Op_Concomitant::init  N = {}", pdata->N);

    target = 0;
    tol0 = 1.0;
    tol = (1.0 - cushion) * tol0;

    spec_norm2 = 1.0;
    spec_norm = 1.0;

    Ax_size = pdata->Naxis * pdata->N;

    if (do_init_weights) {
        obj_weight = 1.0;
        obj_weight *= weight_mod;
    }

    Operator::init();
}

void Op_Concomitant::forward(Eigen::VectorXd &X, Eigen::VectorXd &out) { out = X; }

void Op_Concomitant::transpose(Eigen::VectorXd &X, Eigen::VectorXd &out) { out = X; }

void Op_Concomitant::prox(Eigen::VectorXd &X) {
    spdlog::trace("Starting Op_Concomitant::prox");

    if (do_equil) {
        X.array() /= eq_rows.array();
    }

    double pos = 0.0;
    double neg = 0.0;

    for (int i = 0; i < X.size(); i++) {
        if (pdata->inv_vec(i) > 0) {
            pos += X(i) * X(i);
        } else if (pdata->inv_vec(i) < 0) {
            neg += X(i) * X(i);
        }
    }

    if (pos > neg) {
        double scale = 1.0 / sqrt(pos / neg);
        for (int i = 0; i < X.size(); i++) {
            if (pdata->inv_vec(i) > 0) {
                X(i) *= scale;
            }
        }
    } else if (neg > pos) {
        double scale = 1.0 / sqrt(neg / pos);
        for (int i = 0; i < X.size(); i++) {
            if (pdata->inv_vec(i) < 0) {
                X(i) *= scale;
            }
        }
    }

    if (do_equil) {
        X.array() *= eq_rows.array();
    }

    spdlog::trace("Finished Op_Concomitant::prox");
}

void Op_Concomitant::check(Eigen::VectorXd &X) {
    int is_feas = 1;

    if (do_equil) {
        X.array() /= eq_rows.array();
    }

    double pos = 0.0;
    double neg = 0.0;

    for (int i = 0; i < X.size(); i++) {
        if (pdata->inv_vec(i) > 0) {
            pos += X(i) * X(i);
        } else if (pdata->inv_vec(i) < 0) {
            neg += X(i) * X(i);
        }
    }

    if (abs(1.0 - (pos / neg)) > 0.1) {
        is_feas = 0;
    }

    if (do_equil) {
        X.array() *= eq_rows.array();
    }

    hist_feas.push_back(is_feas);
}

} // namespace Gropt
