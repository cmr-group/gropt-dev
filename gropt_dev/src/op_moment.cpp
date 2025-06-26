#include "spdlog/spdlog.h"

#include "op_moment.hpp"

namespace Gropt {

Op_Moment::Op_Moment(GroptParams &_gparams, double _order)
    : Operator(_gparams)
{
    name = "Moment"; 
    moment_order = _order;
}

Op_Moment::Op_Moment(GroptParams &_gparams, double _order, double _target)
    : Operator(_gparams)
{
    name = "Moment"; 
    moment_order = _order;
    moment_target = _target;
}

void Op_Moment::init()
{
    spdlog::trace("Op_Moment::init  N = {}", gparams->N);
    
    Ax_size = 1;

    A.setZero(1, gparams->Naxis * gparams->N);

    int i_start;
    if (input_start <= 0) {
        i_start = moment_axis*gparams->N;
    } else {
        i_start = input_start + moment_axis*gparams->N;
    }

    int i_stop;
    if (input_stop <= 0) {
        i_stop = (moment_axis + 1)*gparams->N;
    } else {
        i_stop = input_stop + moment_axis*gparams->N;
    }

    spec_norm2 = 0.0;
    for(int j = i_start; j < i_stop; j++) {
        double jj = j - moment_axis*gparams->N;
        double val = 1000.0 * 1000.0 * gparams->dt * pow( (1000.0 * (gparams->dt*(jj - moment_ref0))), moment_order);
        
        A(0, j) = val * gparams->inv_vec(j);
        spec_norm2 += val*val;
    }
    spec_norm2 = sqrt(spec_norm2);
    spec_norm = sqrt(spec_norm2);

    target = moment_target;
    tol0 = moment_tol0;
    tol = (1.0-cushion) * tol0;
    
    if (!save_weights) {
        weight = 1.0e4 / spec_norm2;
    }

    Operator::init();

    spdlog::trace("Initialized operator: {}", name);
    spdlog::trace("    moment_axis = {:d}  moment_order = {:.1f}", moment_axis, moment_order);
    spdlog::trace("    target = {:.1e}  tol0 = {:.1e}  tol = {:.1e}", target, tol0, tol);
    spdlog::trace("    i_start = {:d}  i_stop = {:d}", i_start, i_stop);
}

void Op_Moment::forward(Eigen::VectorXd &X, Eigen::VectorXd &out)
{
    out = A*X;
}

void Op_Moment::transpose(Eigen::VectorXd &X, Eigen::VectorXd &out)
{
    out = A.transpose() * X;
}

void Op_Moment::prox(Eigen::VectorXd &X)
{
    spdlog::trace("Starting Op_Moment::prox");

    for (int i = 0; i < X.size(); i++) {
        double lower_bound = (target-tol);
        double upper_bound = (target+tol);
        X(i) = X(i) < lower_bound ? lower_bound:X(i);
        X(i) = X(i) > upper_bound ? upper_bound:X(i);
    }   

    // for (int i = 0; i < X.size(); i++) {
    //     X(i) = target;
    // }   

    spdlog::trace("Finished Op_Moment::prox");
}


void Op_Moment::check(Eigen::VectorXd &X)
{
    int is_feas = 1;

    for (int i = 0; i < X.size(); i++) {
        double lower_bound = (target-tol0);
        double upper_bound = (target+tol0);

        if ((X(i) < lower_bound) || (X(i) > upper_bound)) {
            is_feas = 0;
        }
    }   

    hist_feas.push_back(is_feas);
}

}  // close "namespace Gropt"