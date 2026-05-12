#include "spdlog/spdlog.h"

#include "op_main.hpp"
#include "workspace_solver.hpp"

namespace Gropt {

Operator::Operator(const ProblemData &_pdata) {
    name = "OperatorMain";
    pdata = &_pdata;
}

Operator::~Operator() {}

void Operator::init() {
    spdlog::trace("In Operator::init() from {}", name);

    N = pdata->N;
    Naxis = pdata->Naxis;
    dt = pdata->dt;
    Ntot = N * Naxis;

    x_temp.setZero(Ntot);
    x_temp_obj.setZero(Ntot);
    Ax_temp.setZero(Ax_size);

    eq_rows.setOnes(Ax_size);
    eq_cols.setOnes(Ntot);
}

void Operator::forward(Eigen::VectorXd &X, Eigen::VectorXd &out) {
    spdlog::warn("Operator::forward is not implemented for the base Operator class.");
}

void Operator::transpose(Eigen::VectorXd &X, Eigen::VectorXd &out) {
    spdlog::warn("Operator::transpose is not implemented for the base Operator class.");
}

void Operator::forward_op(Eigen::VectorXd &X, Eigen::VectorXd &out) {

    if (do_equil) {
        x_temp.array() = X.array() * eq_cols.array();
    } else {
        x_temp = X;
    }

    forward(x_temp, out);

    out.array() /= spec_norm;
    if (do_equil) {
        out.array() *= eq_rows.array();
    }
}

void Operator::transpose_op(Eigen::VectorXd &X, Eigen::VectorXd &out, bool apply_fixer) {
    if (do_equil) {
        Ax_temp.array() = X.array() * eq_rows.array();
    } else {
        Ax_temp = X;
    }

    transpose(Ax_temp, out);

    if (apply_fixer) {
        out.array() *= pdata->fixer.array();
    }
    out.array() /= spec_norm;

    if (do_equil) {
        out.array() *= eq_cols.array();
    }
}

void Operator::transpose_op(Eigen::VectorXd &X, Eigen::VectorXd &out) { transpose_op(X, out, true); }

void Operator::add_Atb(Eigen::VectorXd &b, const WorkspaceSolver &ws) {
    spdlog::trace("Operator::add_Atb  start  name = {}", name);

    Ax_temp.setZero();
    x_temp.setZero();

    Ax_temp = ws.weight * ws.z0 - ws.y0;

    transpose_op(Ax_temp, x_temp);
    b += x_temp;

    spdlog::trace("Operator::add_Atb  end    name = {}", name);
}

void Operator::add_AtAx(Eigen::VectorXd &X, Eigen::VectorXd &out, const WorkspaceSolver &ws) {
    spdlog::trace("Operator::add_AtAx  start  name = {}", name);

    Ax_temp.setZero();
    x_temp.setZero();
    forward_op(X, Ax_temp);
    transpose_op(Ax_temp, x_temp);

    out.array() += ws.weight * x_temp.array();

    spdlog::trace("Operator::add_AtAx  end    name = {}", name);
}

void Operator::add_obj(Eigen::VectorXd &X, Eigen::VectorXd &out) {
    Ax_temp.setZero();
    x_temp.setZero();

    forward_op(X, Ax_temp);
    transpose_op(Ax_temp, x_temp);

    out.array() += obj_weight * x_temp.array();

    spdlog::trace("Operator::add_obj   name = {}  obj_weight = {:.1e}", name, obj_weight);
}

void Operator::check(Eigen::VectorXd &X) {
    feas_check = (X.array() - target).abs().maxCoeff();

    if (feas_check <= tol0) {
        hist_feas.push_back(1);
    } else {
        hist_feas.push_back(0);
    }
}

void Operator::prox(Eigen::VectorXd &X) {
    spdlog::warn("Operator::prox is not implemented for the base Operator class.");
}

void Operator::get_feas(Eigen::VectorXd &s) {
    feas_temp = s;
    prox(feas_temp);
    feas_temp = s - feas_temp;

    r_feas = feas_temp.cwiseAbs().maxCoeff() / (s.cwiseAbs().maxCoeff() + 1.0e-32);

    hist_r_feas.push_back(r_feas);
}

void Operator::print_details() {
    spdlog::info("========== Operator name: {:^16} ==========", name);
    spdlog::info("N: {}, Naxis: {}, Ntot: {}, dt: {}", N, Naxis, Ntot, dt);
    spdlog::info("target: {:.2e}, tol0: {:.2e}, tol: {:.2e}, cushion: {:.2e}", target, tol0, tol, cushion);
    spdlog::info("Ax_size: {}", Ax_size);
    spdlog::info("spec_norm: {:.2e}, spec_norm2: {:.2e}", spec_norm, spec_norm2);
    spdlog::info("obj_weight: {:.2e}, weight_mod: {:.2e}", obj_weight, weight_mod);
    spdlog::info("do_equil: {}, rot_variant: {}, do_init_weights: {}", do_equil, rot_variant, do_init_weights);
    spdlog::default_logger()->flush();
}

} // namespace Gropt
