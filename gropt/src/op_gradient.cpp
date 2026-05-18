#include "spdlog/spdlog.h"

#include "op_gradient.hpp"

namespace Gropt {

Op_Gradient::Op_Gradient(const ProblemData &_pdata, double _gmax, bool _rot_variant, double _weight_mod)
    : Operator(_pdata) {
    name = "Gradient";
    gmax = _gmax;
    rot_variant = _rot_variant;
    weight_mod = _weight_mod;
}

Op_Gradient::Op_Gradient(const ProblemData &_pdata, const Eigen::VectorXd &_gmax_vec, bool _rot_variant,
                         double _weight_mod)
    : Operator(_pdata) {
    if (!_rot_variant) {
        throw std::invalid_argument("Op_Gradient vector gmax is only supported with rot_variant=true");
    }
    name = "Gradient";
    gmax_vec = _gmax_vec;
    gmax = (gmax_vec.size() > 0) ? gmax_vec.maxCoeff() : 0.0;
    rot_variant = _rot_variant;
    weight_mod = _weight_mod;
}

void Op_Gradient::init() {
    spdlog::trace("Op_Gradient::init  N = {}", pdata->N);

    target = 0;
    tol0 = gmax;
    tol = (1.0 - cushion) * tol0;

    spec_norm2 = 1.0;
    spec_norm = 1.0;

    Ax_size = pdata->Naxis * pdata->N;

    if (gmax_vec.size() > 0) {
        int full_size = pdata->Naxis * pdata->N;
        if (gmax_vec.size() == pdata->N && pdata->Naxis > 1) {
            Eigen::VectorXd broadcast(full_size);
            for (int i_ax = 0; i_ax < pdata->Naxis; i_ax++) {
                broadcast.segment(i_ax * pdata->N, pdata->N) = gmax_vec;
            }
            gmax_vec = broadcast;
        } else if (gmax_vec.size() != full_size) {
            throw std::invalid_argument(
                "Op_Gradient: gmax_vec size must be N or Naxis*N");
        }
    }

    if (do_init_weights) {
        obj_weight = 1.0;
        obj_weight *= weight_mod;
    }

    Operator::init();
}

void Op_Gradient::forward(Eigen::VectorXd &X, Eigen::VectorXd &out) { out = X; }

void Op_Gradient::transpose(Eigen::VectorXd &X, Eigen::VectorXd &out) { out = X; }

void Op_Gradient::prox(Eigen::VectorXd &X) {
    spdlog::trace("Starting Op_Gradient::prox");

    if (do_equil) {
        X.array() /= eq_rows.array();
    }

    if (rot_variant) {
        bool use_vec = (gmax_vec.size() == X.size());
        for (int i = 0; i < X.size(); i++) {
            double tol_i = use_vec ? (1.0 - cushion) * gmax_vec(i) : tol;
            double lower_bound = target - tol_i;
            double upper_bound = target + tol_i;
            X(i) = X(i) < lower_bound ? lower_bound : X(i);
            X(i) = X(i) > upper_bound ? upper_bound : X(i);

            // This is specific to the Op_Gradient operator
            if (!isnan(pdata->set_vals(i))) {
                X(i) = pdata->set_vals(i);
            }
        }
    } else {
        for (int i = 0; i < N; i++) {
            double upper_bound = (target + tol);

            double val = 0.0;
            for (int i_ax = 0; i_ax < Naxis; i_ax++) {
                val += X(i_ax * N + i) * X(i_ax * N + i);
            }
            val = sqrt(val);

            if (val > upper_bound) {
                for (int i_ax = 0; i_ax < Naxis; i_ax++) {
                    X(i_ax * N + i) *= (upper_bound / val);
                }
            }
        }

        for (int i = 0; i < X.size(); i++) {
            if (!isnan(pdata->set_vals(i))) {
                X(i) = pdata->set_vals(i);
            }
        }
    }

    if (do_equil) {
        X.array() *= eq_rows.array();
    }

    spdlog::trace("Finished Op_Gradient::prox");
}

void Op_Gradient::check(Eigen::VectorXd &X) {
    int is_feas = 1;

    if (do_equil) {
        X.array() /= eq_rows.array();
    }

    if (rot_variant) {
        bool use_vec = (gmax_vec.size() == X.size());
        for (int i = 0; i < X.size(); i++) {
            double tol_i = use_vec ? gmax_vec(i) : tol0;
            double lower_bound = target - tol_i;
            double upper_bound = target + tol_i;

            if ((X(i) < lower_bound) || (X(i) > upper_bound) && isnan(pdata->set_vals(i))) {
                is_feas = 0;
            }
        }
    } else {
        for (int i = 0; i < N; i++) {
            double upper_bound = (target + tol0);

            double val = 0.0;
            for (int i_ax = 0; i_ax < Naxis; i_ax++) {
                val += X(i_ax * N + i) * X(i_ax * N + i);
            }
            val = sqrt(val);

            if ((val > upper_bound) && isnan(pdata->set_vals(i))) {
                is_feas = 0;
            }
        }
    }

    if (do_equil) {
        X.array() *= eq_rows.array();
    }

    hist_feas.push_back(is_feas);
}

} // namespace Gropt
