#include "spdlog/spdlog.h"

#include "op_moment.hpp"

namespace Gropt {

Op_Moment::Op_Moment(GroptParams &_gparams, double _order, double _target, double _tol0, std::string _units,
                    int _moment_axis, int _start_idx0, int _stop_idx0, int _ref_idx0, double _weight_mod)
    : Operator(_gparams)
{
    name = "Moment"; 
    moment_order = _order;
    units = _units;

    double moment_scale = 1.0;
    if (units == "mT*ms/m") {
        moment_scale = 1.0;  
    } else if (units == "T*s/m") {
        moment_scale = 1000.0 * pow(1000.0, moment_order+1);   
    } else if (units == "rad*s/m") {
        moment_scale = 1000.0 * pow(1000.0, moment_order+1) / 4.257638544e7; 
    } else if (units == "s/m") {
        moment_scale = 1000.0 * pow(1000.0, moment_order+1) / 2.675153194e8; 
    } else {
        spdlog::error("Unsupported units for moment constraint: {}", units);
        throw std::invalid_argument("Unsupported units for moment constraint");
    }

    moment_target = _target * moment_scale;
    moment_tol0 = _tol0 * moment_scale;
    moment_axis = _moment_axis;


    start_idx0 = _start_idx0;
    stop_idx0 = _stop_idx0;
    ref_idx0 = _ref_idx0;   

    start_idx = start_idx0;
    stop_idx = stop_idx0;
    ref_idx = ref_idx0;

    weight_mod = _weight_mod;
}

void Op_Moment::init()
{
    spdlog::trace("Op_Moment::init  N = {}", gparams->N);
    
    Ax_size = 1;

    A.setZero(1, gparams->Naxis * gparams->N);

    // If start and stop indices are not set, constraint covers the whole axis
    int i_start;
    if (start_idx <= 0) {
        i_start = moment_axis*gparams->N;
    } else {
        i_start = start_idx + moment_axis*gparams->N;
    }

    int i_stop;
    if (stop_idx <= 0) {
        i_stop = (moment_axis + 1)*gparams->N;
    } else {
        i_stop = stop_idx + moment_axis*gparams->N;
    }

    spec_norm2 = 0.0;
    for(int j = i_start; j < i_stop; j++) {
        double jj = j - moment_axis*gparams->N;
        double val = 1000.0 * 1000.0 * gparams->dt * pow( (1000.0 * (gparams->dt*(jj - ref_idx))), moment_order);
        
        A(0, j) = val * gparams->inv_vec(j);
        spec_norm2 += val*val;
    }
    // TODO: I think this sqrt is wrong, only the second one is needed, do some tests to confirm
    spec_norm2 = sqrt(spec_norm2);
    spec_norm = sqrt(spec_norm2);

    target = moment_target;
    tol0 = moment_tol0;
    tol = (1.0-cushion) * tol0;
    
    if (do_init_weights) {
        weight = 1.0e4;
        obj_weight = 1.0;

        weight *= weight_mod;
        obj_weight *= weight_mod;
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

    // It can work to set values to the target and not tolerance, TODO: add an option for this
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