#include "spdlog/spdlog.h"

#include "op_tv.hpp"

namespace Gropt {

Op_TV::Op_TV(const ProblemData &_pdata, double _tv_lam, int _order, double _weight_mod)
    : Operator(_pdata)
{
    name = "TotalVariation";
    weight_mod = _weight_mod;
    tv_lam = _tv_lam;
    order = _order;
}

void Op_TV::init()
{
    spdlog::trace("Op_TV::init  N = {}", pdata->N);

    target = 0;
    tol0 = tv_lam;
    tol = (1.0-cushion) * tol0;

    spec_norm2 = std::pow(4.0/(pdata->dt*pdata->dt), order);
    spec_norm = sqrt(spec_norm2);

    if (order < 0 || order > pdata->N) {
        spdlog::error("Op_TV initialized with invalid order {}!", order);
    }
    Ax_size = pdata->Naxis * (pdata->N-order);

    if (do_init_weights) {
        obj_weight = 1.0;
        obj_weight *= weight_mod;
    }

    Operator::init();
}

void Op_TV::forward(Eigen::VectorXd &X, Eigen::VectorXd &out)
{
    if (order == 0) {
        out = X;
        return;
    }

    Eigen::VectorXd temp_in = X;
    Eigen::VectorXd temp_out(Naxis * (pdata->N-1));
    int curr_len = pdata->N;

    for (int it = 0; it < order; it++) {
        temp_out.setZero(Naxis * (curr_len - 1));
        for (int i_ax = 0; i_ax < Naxis; i_ax++) {
            for (int i = 0; i < (curr_len - 1); i++) {
                temp_out(i_ax*(curr_len - 1)+i) = (temp_in(i_ax*curr_len+i+1) - temp_in(i_ax*curr_len+i)) / pdata->dt;
            }
        }
        curr_len--;
        temp_in = temp_out;
    }
    out = temp_out;
}

void Op_TV::transpose(Eigen::VectorXd &X, Eigen::VectorXd &out)
{
    if (order == 0) {
        out = X;
        return;
    }

    Eigen::VectorXd temp_in = X;
    Eigen::VectorXd temp_out(Naxis * (pdata->N - order + 1));
    int curr_len = pdata->N - order;

    for (int it = order; it > 0; it--) {
        temp_out.setZero(Naxis * (curr_len + 1));
        for (int i_ax = 0; i_ax < Naxis; i_ax++) {
            temp_out(i_ax*(curr_len+1)+0) = -temp_in(i_ax*curr_len+0) / pdata->dt;
            for (int i = 1; i < curr_len; i++) {
                temp_out(i_ax*(curr_len+1)+i) = (temp_in(i_ax*curr_len+i-1) - temp_in(i_ax*curr_len+i)) / pdata->dt;
            }
            temp_out(i_ax*(curr_len+1)+curr_len) = temp_in(i_ax*curr_len+curr_len-1) / pdata->dt;
        }
        curr_len++;
        temp_in = temp_out;
    }
    out = temp_out;
}

void Op_TV::prox(Eigen::VectorXd &X)
{
    spdlog::trace("Starting Op_TV::prox");

    for (int i = 0; i < X.size(); i++) {
        if (abs(X(i)) > tv_lam) {
            X(i) = X(i) > 0 ? X(i) - tv_lam : X(i) + tv_lam;
        } else {
            X(i) = 0.0;
        }
    }

    spdlog::trace("Finished Op_TV::prox");
}


void Op_TV::check(Eigen::VectorXd &X)
{
    int is_feas = 1;

    // As of right now we will assume that the TV operator is always feasible
    // This isn't necessarily the case, we should have an option to check feasibility (TODO)

    hist_feas.push_back(is_feas);
}

}  // close "namespace Gropt"
