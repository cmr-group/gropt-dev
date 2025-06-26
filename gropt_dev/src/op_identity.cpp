#include "spdlog/spdlog.h"

#include "op_identity.hpp"

namespace Gropt {

Op_Identity::Op_Identity(GroptParams &_gparams)
    : Operator(_gparams)
{
    name = "Identity"; 
}

Op_Identity::Op_Identity(GroptParams &_gparams, double _weight_mod)
    : Operator(_gparams)
{
    name = "Identity"; 
    weight_mod = _weight_mod;
}

void Op_Identity::init()
{
    spdlog::trace("Op_Identity::init  N = {}", gparams->N);
    
    target = 0;

    spec_norm2 = 1.0;
    spec_norm = 1.0;

    Ax_size = gparams->Naxis * gparams->N;

    weight = weight_mod;
    obj_weight = weight_mod;

    Operator::init();
}

void Op_Identity::forward(Eigen::VectorXd &X, Eigen::VectorXd &out)
{
    out = X;
}

void Op_Identity::transpose(Eigen::VectorXd &X, Eigen::VectorXd &out)
{
    out = X;    
}

void Op_Identity::prox(Eigen::VectorXd &X)
{
    spdlog::trace("Starting Op_Identity::prox");
    spdlog::trace("Finished Op_Identity::prox");
}


void Op_Identity::check(Eigen::VectorXd &X)
{
    int is_feas = 1;
    hist_feas.push_back(is_feas);
}

}  // close "namespace Gropt"