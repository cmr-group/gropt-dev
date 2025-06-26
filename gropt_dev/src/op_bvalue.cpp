#include "spdlog/spdlog.h"

#include "op_bvalue.hpp"

namespace Gropt {

Op_BValue::Op_BValue(GroptParams &_gparams)
    : Operator(_gparams)
{
    name = "b-value"; 
}

void Op_BValue::init()
{
    spdlog::trace("Op_BValue::init  N = {}", gparams->N);
    
    tol = (1.0-cushion) * tol0;

    GAMMA = 267.5221900e6;  // rad/S/T
    MAT_SCALE = pow((GAMMA / 1000.0 * gparams->dt), 2.0) * gparams->dt;  // 1/1000 is for m->mm in b-value

    spec_norm2 = (gparams->N*gparams->N+gparams->N)/2.0 * MAT_SCALE;;
    spec_norm = sqrt(spec_norm2);

    Ax_size = gparams->Naxis * gparams->N;

    if (!save_weights) {
        weight = 1.0e4 / spec_norm2;
        obj_weight = -1.0/spec_norm2;
    }

    Operator::init();
}

void Op_BValue::forward(Eigen::VectorXd &X, Eigen::VectorXd &out)
{
    out.setZero();
    for (int j = 0; j < Naxis; j++) {
        int jN = j*N;
        double gt = 0;    
        for (int i = start_ind; i < N; i++) {
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
        for (int i = N-1; i >= start_ind; i--) {
            gt += X(jN + i) * sqrt(MAT_SCALE);
            out(jN + i) = gt * gparams->inv_vec(jN + i);
        }
    }
}

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