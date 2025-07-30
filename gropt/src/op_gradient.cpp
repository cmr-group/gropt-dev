#include "spdlog/spdlog.h"

#include "op_gradient.hpp"

namespace Gropt {

Op_Gradient::Op_Gradient(GroptParams &_gparams, double _gmax, bool _rot_variant, double _weight_mod)
    : Operator(_gparams)
{
    name = "Gradient"; 
    gmax = _gmax;
    rot_variant = _rot_variant;
    weight_mod = _weight_mod;
}

void Op_Gradient::init()
{
    spdlog::trace("Op_Gradient::init  N = {}", gparams->N);
    
    target = 0;
    tol0 = gmax;
    tol = (1.0-cushion) * tol0;

    spec_norm2 = 1.0;
    spec_norm = 1.0;

    Ax_size = gparams->Naxis * gparams->N;

    if (do_init_weights) {
        weight = 1.0;
        obj_weight = 1.0;

        weight *= weight_mod;
        obj_weight *= weight_mod;
    }

    Operator::init();
}

void Op_Gradient::forward(Eigen::VectorXd &X, Eigen::VectorXd &out)
{
    out = X;
}

void Op_Gradient::transpose(Eigen::VectorXd &X, Eigen::VectorXd &out)
{
    out = X;    
}

void Op_Gradient::prox(Eigen::VectorXd &X)
{
    spdlog::trace("Starting Op_Gradient::prox");

    if (rot_variant) {
        for (int i = 0; i < X.size(); i++) {
            double lower_bound = (target-tol);
            double upper_bound = (target+tol);
            X(i) = X(i) < lower_bound ? lower_bound:X(i);
            X(i) = X(i) > upper_bound ? upper_bound:X(i);
            
            // This is specific to the Op_Gradient operator
            if (!isnan(gparams->set_vals(i))) {
                X(i) = gparams->set_vals(i);
            }
        }   
    } else {
        for (int i = 0; i < N; i++) {
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

        for (int i = 0; i < X.size(); i++) {
            if (!isnan(gparams->set_vals(i))) {
                X(i) = gparams->set_vals(i);
            }
        }
    }

    spdlog::trace("Finished Op_Gradient::prox");
}


void Op_Gradient::check(Eigen::VectorXd &X)
{
    int is_feas = 1;

    if (rot_variant) {
        for (int i = 0; i < X.size(); i++) {
            double lower_bound = (target-tol0);
            double upper_bound = (target+tol0);

            if ((X(i) < lower_bound) || (X(i) > upper_bound) && isnan(gparams->set_vals(i))) {
                is_feas = 0;
            }
        }   
    } else {
        for (int i = 0; i < N; i++) {
            double upper_bound = (target+tol0);
            
            double val = 0.0;
            for (int i_ax = 0; i_ax < Naxis; i_ax++) {
                val += X(i_ax*N+i)*X(i_ax*N+i);
            }
            val = sqrt(val);

            if ((val > upper_bound) && isnan(gparams->set_vals(i))) {
                is_feas = 0;
            }
        }
    }

    hist_feas.push_back(is_feas);

}

}  // close "namespace Gropt"