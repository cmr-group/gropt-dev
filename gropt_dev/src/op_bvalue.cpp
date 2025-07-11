#include "spdlog/spdlog.h"

#include "op_bvalue.hpp"

namespace Gropt {

Op_BValue::Op_BValue(GroptParams &_gparams, double _bval_target, double _bval_tol0,
                     int _start_idx0, int _stop_idx0, double _weight_mod)
    : Operator(_gparams)
{
    name = "b-value"; 

    bval_target = _bval_target;
    bval_tol0 = _bval_tol0;

    start_idx0 = _start_idx0;
    stop_idx0 = _stop_idx0;

    start_idx = start_idx0;
    stop_idx = stop_idx0;

    weight_mod = _weight_mod;
}

void Op_BValue::init()
{
    spdlog::trace("Op_BValue::init  N = {}", gparams->N);
    
    target = bval_target;
    tol0 = bval_tol0;
    tol = (1.0-cushion) * tol0;

    GAMMA = 267.5221900e6;  // rad/S/T
    MAT_SCALE = pow((GAMMA / 1000.0 * gparams->dt), 2.0) * gparams->dt;  // 1/1000 is for m->mm in b-value

    // If start and stop indices are not set, constraint covers the whole axis
    if (start_idx <= 0) {
        i_start = 0;
    } else {
        i_start = start_idx;
    }

    if (stop_idx <= 0) {
        i_stop = gparams->N;
    } else {
        i_stop = stop_idx;
    }

    int Nnorm = i_stop - i_start;
    spec_norm2 = (Nnorm*Nnorm + Nnorm)/2.0 * MAT_SCALE;
    spec_norm = sqrt(spec_norm2);

    Ax_size = gparams->Naxis * gparams->N;

    if (do_init_weights) {
        weight = 1.0e4;
        obj_weight = -1.0;

        weight *= weight_mod;
        obj_weight *= weight_mod;
    }

    Operator::init();
}

void Op_BValue::forward(Eigen::VectorXd &X, Eigen::VectorXd &out)
{
    out.setZero();
    for (int j = 0; j < Naxis; j++) {
        int jN = j*N;
        double gt = 0;    
        for (int i = i_start; i < i_stop; i++) {
            gt += X(jN + i) * gparams->inv_vec(jN + i);
            out(jN + i) = gt * sqrt(MAT_SCALE);
        }
    }
}

void Op_BValue::transpose(Eigen::VectorXd &X, Eigen::VectorXd &out)
{
    out.setZero(); 
    for (int j = 0; j < Naxis; j++) {
        int jN = j*N;
        double gt = 0;    
        for (int i = i_stop-1; i >= i_start; i--) {
            gt += X(jN + i) * sqrt(MAT_SCALE);
            out(jN + i) = gt * gparams->inv_vec(jN + i);
        }
    }
}

// TODO: This seems all wrong for three axis case, need to fix
void Op_BValue::prox(Eigen::VectorXd &X)
{
    spdlog::trace("Starting Op_BValue::prox");
    for (int j = 0; j < Naxis; j++) {
        double xnorm = X.segment(j*N, N).norm();
        double min_val = sqrt(target - tol);
        double max_val = sqrt(target + tol);

        if (xnorm < min_val) {
            X.segment(j*N, N) *= (min_val/xnorm);
        } else if (xnorm > max_val) {
            X.segment(j*N, N) *= (max_val/xnorm);
        }
    }

    spdlog::trace("Finished Op_BValue::prox");
}


void Op_BValue::check(Eigen::VectorXd &X)
{
    int is_feas = 1;

    for (int j = 0; j < Naxis; j++) {
        double bval_t = (X.segment(j*N, N)).squaredNorm();    
        
        double d_bval = fabs(bval_t - target);
        if (d_bval > tol0) {
            is_feas = 0;
        }

    }

    hist_feas.push_back(is_feas);
}

double Op_BValue::get_bvalue(Eigen::VectorXd &X)
{
    Ax_temp.setZero();
    forward_op(X, Ax_temp);

    return Ax_temp.squaredNorm();
}

}  // close "namespace Gropt"