#include "spdlog/spdlog.h"

#include "op_tv.hpp"

namespace Gropt {

Op_TV::Op_TV(GroptParams &_gparams, double _tv_lam, double _weight_in)
    : Operator(_gparams)
{
    name = "TotalVariation"; 
    weight_in = _weight_in;
    tv_lam = _tv_lam;
}

void Op_TV::init()
{
    spdlog::trace("Op_TV::init  N = {}", gparams->N);
    
    // These arent really used for shrinkage
    target = 0;
    tol0 = tv_lam;
    tol = (1.0-cushion) * tol0;

    spec_norm2 = 4.0/gparams->dt/gparams->dt;
    spec_norm = sqrt(spec_norm2);

    Ax_size = gparams->Naxis * (gparams->N-1);

    if (!save_weights) {
        weight = weight_in * 1e4;
    }

    Operator::init();
}

void Op_TV::forward(Eigen::VectorXd &X, Eigen::VectorXd &out)
{
    for (int i_ax = 0; i_ax < Naxis; i_ax++) {
        for (int i = 0; i < (N-1); i++) {
            out(i_ax*(N-1)+i) = (X(i_ax*N+i+1) - X(i_ax*N+i))/gparams->dt;
        }
    }
}

void Op_TV::transpose(Eigen::VectorXd &X, Eigen::VectorXd &out)
{
    for (int i_ax = 0; i_ax < Naxis; i_ax++) {
        out(i_ax*N+0) = -X(i_ax*(N-1)+0) / gparams->dt;
        for (int i = 1; i < (N-1); i++) {
            out(i_ax*N+i) = (X(i_ax*(N-1)+i-1) - X(i_ax*(N-1)+i)) / gparams->dt;
        }
        out(i_ax*N+N-1) = X(i_ax*(N-1)+N-2) / gparams->dt;
    }
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