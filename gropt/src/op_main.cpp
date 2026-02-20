#include "spdlog/spdlog.h"

#include "op_main.hpp"
#include "sdmm_workspace.hpp"

namespace Gropt {

    Operator::Operator(const ProblemData &_pdata)
    {
        name = "OperatorMain";
        pdata = &_pdata;
    }

    Operator::~Operator() {}

    void Operator::init()
    {
        spdlog::trace("In Operator::init() from {}", name);

        N = pdata->N;
        Naxis = pdata->Naxis;
        dt = pdata->dt;
        Ntot = N * Naxis;

        x_temp.setZero(Ntot);
        x_temp_obj.setZero(Ntot);
        Ax_temp.setZero(Ax_size);
    }

    void Operator::forward(Eigen::VectorXd &X, Eigen::VectorXd &out)
    {
        spdlog::warn("Operator::forward is not implemented for the base Operator class.");
    }

    void Operator::transpose(Eigen::VectorXd &X, Eigen::VectorXd &out)
    {
        spdlog::warn("Operator::transpose is not implemented for the base Operator class.");
    }

    void Operator::forward_op(Eigen::VectorXd &X, Eigen::VectorXd &out)
    {
        forward(X, out);
    }

    void Operator::transpose_op(Eigen::VectorXd &X, Eigen::VectorXd &out)
    {
        transpose(X, out);
        out.array() *= pdata->fixer.array();
        out.array() /= spec_norm2;
    }

    void Operator::add_Atb(Eigen::VectorXd &b, const SDMMWorkspace &ws)
    {
        spdlog::trace("Operator::add_Atb  start  name = {}", name);

        Ax_temp.setZero();
        x_temp.setZero();

        Ax_temp = ws.weight*ws.z0 - ws.y0;

        transpose_op(Ax_temp, x_temp);
        b += x_temp;

        spdlog::trace("Operator::add_Atb  end    name = {}", name);
    }

    void Operator::add_AtAx(Eigen::VectorXd &X, Eigen::VectorXd &out, const SDMMWorkspace &ws)
    {
        spdlog::trace("Operator::add_AtAx  start  name = {}", name);

        Ax_temp.setZero();
        x_temp.setZero();
        forward_op(X, Ax_temp);
        transpose_op(Ax_temp, x_temp);

        out.array() += ws.weight*x_temp.array();

        spdlog::trace("Operator::add_AtAx  end    name = {}", name);
    }

    void Operator::add_obj(Eigen::VectorXd &X, Eigen::VectorXd &out)
    {
        Ax_temp.setZero();
        x_temp.setZero();

        forward_op(X, Ax_temp);
        transpose_op(Ax_temp, x_temp);

        out.array() += obj_weight*Ax_temp.array();

        spdlog::trace("Operator::add_obj   name = {}  obj_weight = {:.1e}", name, obj_weight);
    }

    void Operator::check(Eigen::VectorXd &X)
    {
        feas_check = (X.array() - target).abs().maxCoeff();

        if (feas_check <= tol0) {
            hist_feas.push_back(1);
        } else {
            hist_feas.push_back(0);
        }
    }

    void Operator::prox(Eigen::VectorXd &X)
    {
        spdlog::warn("Operator::prox is not implemented for the base Operator class.");
    }

    void Operator::get_feas(Eigen::VectorXd &s)
    {
        feas_temp = s;
        prox(feas_temp);
        feas_temp = s - feas_temp;

        r_feas = feas_temp.cwiseAbs().maxCoeff()/(s.cwiseAbs().maxCoeff() + 1.0e-32);

        hist_r_feas.push_back(r_feas);
    }

}  // close "namespace Gropt"
