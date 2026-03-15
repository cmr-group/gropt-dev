#include "spdlog/spdlog.h"

#include "op_eddy.hpp"

// TODO:  When we build A with Naxis > 1, A could just have full length so forward and transpose are
//        simple matrix multiplications instead of the segmented thing now.

namespace Gropt {
Op_Eddy::Op_Eddy(const ProblemData &_pdata, const Eigen::VectorXd &_lam, double _tol, double _weight_mod)
    : Operator(_pdata) {
    name = "Eddy";

    weight_mod = _weight_mod;
    tol0 = _tol;
    lambdas = _lam;

    Nlam = lambdas.size();
}

void Op_Eddy::init() {

    Nrows = Nlam * pdata->Naxis;
    target = 0;
    tol = (1.0 - cushion) * tol0;

    A.setZero(Nrows, pdata->N);

    spec_norm2 = 0.0;
    for (int i = 0; i < pdata->Naxis; i++) {
        for (int i_lam = 0; i_lam < Nlam; i_lam++) {

            int ii = i_lam + i * Nlam;
            double lam = lambdas(i_lam);

            for (int j = 0; j < pdata->N; j++) {
                double jj = pdata->N - j - 1;
                double val = exp(-(jj + 1.0) * pdata->dt / lam) - exp(-jj * pdata->dt / lam);
                A(ii, j) = val;
            }

            spec_norm2 += A.row(ii).squaredNorm();
        }
    }

    spec_norm = sqrt(spec_norm2);
    Ax_size = Nrows;

    Operator::init();
}

void Op_Eddy::forward(Eigen::VectorXd &X, Eigen::VectorXd &out) {
    // TODO: Fix the Naxis Nlam selection
    for (int i = 0; i < Naxis; i++) {
        for (int i_lam = 0; i_lam < Nlam; i_lam++) {
            int ii = i_lam + i * Nlam;
            out(ii) = A.row(ii) * X.segment(i * N, N);
        }
    }
}

void Op_Eddy::transpose(Eigen::VectorXd &X, Eigen::VectorXd &out) {
    out.setZero();

    for (int i = 0; i < Naxis; i++) {
        for (int i_lam = 0; i_lam < Nlam; i_lam++) {
            int ii = i_lam + i * Nlam;
            out.segment(i * N, N) += A.row(ii).transpose() * X(ii);
        }
    }
}

void Op_Eddy::prox(Eigen::VectorXd &X) {
    spdlog::trace("Starting Op_Eddy::prox");

    if (do_equil) {
        X.array() /= eq_rows.array();
    }
    X.array() *= spec_norm;

    for (int i = 0; i < X.size(); i++) {
        double lower_bound = (target - tol);
        double upper_bound = (target + tol);
        X(i) = X(i) < lower_bound ? lower_bound : X(i);
        X(i) = X(i) > upper_bound ? upper_bound : X(i);
    }

    if (do_equil) {
        X.array() *= eq_rows.array();
    }
    X.array() /= spec_norm;
}

void Op_Eddy::check(Eigen::VectorXd &X) {
    int is_feas = 1;

    if (do_equil) {
        X.array() /= eq_rows.array();
    }
    X.array() *= spec_norm;

    for (int i = 0; i < X.size(); i++) {
        double lower_bound = (target - tol0);
        double upper_bound = (target + tol0);

        if ((X(i) < lower_bound) || (X(i) > upper_bound)) {
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