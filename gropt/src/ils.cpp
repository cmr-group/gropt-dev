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

void IndirectLinearSolver::get_lhs(Eigen::VectorXd &x, Eigen::VectorXd &out)
{
    out.setZero();

    if (tik_lam > 0.0) {
        out.array() += tik_lam * x.array();
    }

    out.array() += sigma * x.array();

    for (int i = 0; i < gparams->all_op.size(); i++) {
        if (gparams->all_op[i]->use_projection) continue; // Ignore constraints in project mode
        gparams->all_op[i]->add_AtAx(x, out, *ws[i]);
    }

    // Every objective is offered the LHS; each decides what it adds: convex objectives add AᵀA
    // curvature, DCA/linearized ones add nothing here (they go to the RHS), and augmented-Lagrangian
    // constraint-objectives add their (rank-1 PSD) Gauss-Newton curvature.
    for (int i = 0; i < gparams->all_obj.size(); i++) {
        gparams->all_obj[i]->add_obj(x, out);
    }
}

void IndirectLinearSolver::get_rhs(Eigen::VectorXd &x0, Eigen::VectorXd &out)
{
    out.setZero();

    out.array() += sigma * x0.array();

    spdlog::trace("IndirectLinearSolver::get_rhs  iteration size = {:d}", gparams->all_op.size());
    for (int i = 0; i < gparams->all_op.size(); i++) {
        if (gparams->all_op[i]->use_projection) continue; // Ignore constraints in project mode
        gparams->all_op[i]->add_Atb(out, *ws[i]);
    }

    // Linearized (DCA) objective contribution, frozen at the current iterate x0. Per-objective:
    // only objectives flagged linearize_obj go here; convex ones stay as curvature in get_lhs.
    // normalize_obj self-normalizes the objective direction so weight_mod sets a CONSTANT step
    // magnitude (the pull no longer grows with ||AᵀA x||) — a forgiving rate knob.
    // Every objective is offered the RHS; each decides what it adds: DCA/linearized objectives add
    // their frozen-gradient pull, convex ones add nothing here (they went to the LHS), and
    // augmented-Lagrangian constraint-objectives add their (scalar-dual) pull.
    for (int i = 0; i < gparams->all_obj.size(); i++) {
        gparams->all_obj[i]->add_obj_rhs(x0, out, gparams->normalize_obj);
    }
}

} // close "namespace Gropt"
