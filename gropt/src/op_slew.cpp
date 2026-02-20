#include "spdlog/spdlog.h"

#include "op_slew.hpp"

namespace Gropt {

Op_Slew::Op_Slew(const ProblemData &_pdata, double _smax, bool _rot_variant, double _weight_mod)
    : Operator(_pdata)
{
    name = "Slew";
    smax = _smax;
    rot_variant = _rot_variant;
    weight_mod = _weight_mod;
}

void Op_Slew::init()
{
    spdlog::trace("Op_Slew::init  N = {}", pdata->N);

    target = 0;
    tol0 = smax;
    tol = (1.0-cushion) * tol0;

    spec_norm2 = 4.0/pdata->dt/pdata->dt;
    spec_norm = sqrt(spec_norm2);

    Ax_size = pdata->Naxis * (pdata->N-1);

    if (do_init_weights) {
        obj_weight = 1.0;
        obj_weight *= weight_mod;
    }

    Operator::init();
}

void Op_Slew::forward(Eigen::VectorXd &X, Eigen::VectorXd &out)
{
    for (int i_ax = 0; i_ax < Naxis; i_ax++) {
        for (int i = 0; i < (N-1); i++) {
            out(i_ax*(N-1)+i) = (X(i_ax*N+i+1) - X(i_ax*N+i))/pdata->dt;
        }
    }
}

void Op_Slew::transpose(Eigen::VectorXd &X, Eigen::VectorXd &out)
{
    for (int i_ax = 0; i_ax < Naxis; i_ax++) {
        out(i_ax*N+0) = -X(i_ax*(N-1)+0) / pdata->dt;
        for (int i = 1; i < (N-1); i++) {
            out(i_ax*N+i) = (X(i_ax*(N-1)+i-1) - X(i_ax*(N-1)+i)) / pdata->dt;
        }
        out(i_ax*N+N-1) = X(i_ax*(N-1)+N-2) / pdata->dt;
    }
}

void Op_Slew::prox(Eigen::VectorXd &X)
{
    spdlog::trace("Starting Op_Slew::prox");

    if (rot_variant) {
        for (int i = 0; i < X.size(); i++) {
            double lower_bound = (target-tol);
            double upper_bound = (target+tol);
            X(i) = X(i) < lower_bound ? lower_bound:X(i);
            X(i) = X(i) > upper_bound ? upper_bound:X(i);

        }
    } else {
        for (int i = 0; i < N-1; i++) {
            double upper_bound = (target+tol);

            double val = 0.0;
            for (int i_ax = 0; i_ax < Naxis; i_ax++) {
                val += X(i_ax*N+i)*X(i_ax*N+i);
            }
            val = sqrt(val);

            if (val > upper_bound) {
                for (int i_ax = 0; i_ax < Naxis; i_ax++) {
                    X(i_ax*N+i) *= (upper_bound/val);
                }
            }
        }
    }

    spdlog::trace("Finished Op_Slew::prox");
}


void Op_Slew::check(Eigen::VectorXd &X)
{
    int is_feas = 1;

    if (rot_variant) {
        for (int i = 0; i < X.size(); i++) {

            bool should_check = isnan(pdata->set_vals(i));
            if (i > 0) {should_check = should_check || isnan(pdata->set_vals(i-1)); }
            if (i < X.size()-1) {should_check = should_check || isnan(pdata->set_vals(i+1)); }
            double lower_bound = (target-tol0);
            double upper_bound = (target+tol0);

            if (((X(i) < lower_bound) || (X(i) > upper_bound)) && should_check) {
                is_feas = 0;
            }
        }
    } else {
        for (int i = 0; i < N-1; i++) {
            bool should_check = isnan(pdata->set_vals(i));
            if (i > 0) {should_check = should_check || isnan(pdata->set_vals(i-1)); }
            if (i < N-1) {should_check = should_check || isnan(pdata->set_vals(i+1)); }

            double upper_bound = (target+tol0);

            double val = 0.0;
            for (int i_ax = 0; i_ax < Naxis; i_ax++) {
                val += X(i_ax*N+i)*X(i_ax*N+i);
            }
            val = sqrt(val);

            if ((val > upper_bound) && should_check) {
                is_feas = 0;
            }
        }
    }

    hist_feas.push_back(is_feas);
}

}  // close "namespace Gropt"
