#include "spdlog/spdlog.h"

#include "ils.hpp"

namespace Gropt {

IndirectLinearSolver::IndirectLinearSolver(GroptParams &_gparams, int _n_iter, double _sigma, double _tik_lam)
    : name("IndirectLinearSolver"), 
    n_iter(_n_iter), 
    sigma(_sigma), 
    tik_lam(_tik_lam)
{
    gparams = &_gparams;
    hist_n_iter.push_back(-1);
} 

Eigen::VectorXd IndirectLinearSolver::solve(Eigen::VectorXd &x0)
{
    spdlog::warn("IndirectLinearSolver::solve is not implemented for the base class.");
    return Eigen::VectorXd::Zero(1);
}

Eigen::VectorXd IndirectLinearSolver::get_lhs(Eigen::VectorXd &x, Eigen::VectorXd &out)
{
    out.setZero();

    if (tik_lam > 0.0) {
        out.array() += tik_lam * x.array();
    }

    out.array() += sigma * x.array();

    for (int i = 0; i < gparams->all_op.size(); i++) {
        gparams->all_op[i]->add_AtAx(x, out);
    }

    for (int i = 0; i < gparams->all_obj.size(); i++) {
        gparams->all_obj[i]->add_obj(x, out);
    }

    return out; 
} 

Eigen::VectorXd IndirectLinearSolver::get_rhs(Eigen::VectorXd &x0, Eigen::VectorXd &out)
{
    out.setZero();

    out.array() += sigma * x0.array();

    spdlog::trace("IndirectLinearSolver::get_rhs  iteration size = {:d}", gparams->all_op.size());
    for (int i = 0; i < gparams->all_op.size(); i++) {
        gparams->all_op[i]->add_Atb(out);
    }

    // for (int i = 0; i < gparams->all_obj.size(); i++) {
    //     gparams->all_obj[i]->add_b(out);
    // }

    return out;
}

} // close "namespace Gropt"