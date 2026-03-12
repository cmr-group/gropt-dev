#include "spdlog/spdlog.h"

#include "ils.hpp"
#include "ils_bicgstabl.hpp"
#include "ils_cg.hpp"
#include "ils_nlcg.hpp"
#include "solver_osqp.hpp"

namespace Gropt {

SolveResult SolverOSQP::solve(GroptParams &_gparams) {
    spdlog::trace("Starting SolverOSQP::solve");
    gparams = &_gparams;
    if (gparams->op_prep_status != gparams->N) {
        spdlog::info("Operators do not seem prepared, calling prepare()");
        gparams->prepare();
    }

    // Initialize per-operator OSQP workspaces
    osqp_ws.resize(gparams->all_op.size());
    for (int i = 0; i < gparams->all_op.size(); i++) {
        Operator *op = gparams->all_op[i].get();

        // Set initial weight based on operator type
        osqp_ws[i].weight = 1.0;
        // Slew, moment, bvalue, SAFE, TV operators start with higher weight
        if (op->name == "Slew" || op->name == "Moment" || op->name == "b-value" || op->name == "SAFE" ||
            op->name == "TotalVariation") {
            osqp_ws[i].weight = 1e4;
        }
        osqp_ws[i].weight *= op->weight_mod;

        osqp_ws[i].init(op->Ax_size);
        osqp_ws[i].prep(*op, gparams->pdata.X0);
    }

    // Populate base class ws pointers
    ws.resize(osqp_ws.size());
    for (int i = 0; i < osqp_ws.size(); i++) {
        ws[i] = &osqp_ws[i];
    }

    Eigen::VectorXd X = gparams->pdata.X0;
    if (gparams->all_op[0]->do_equil) {
        X.array() /= gparams->all_op[0]->eq_cols.array();
    }
    Eigen::VectorXd Xhat;

    Px.setZero(X.size());
    r_dual.setZero(X.size());

    if (gparams->ils_method == CG) {
        ils_solver = new ILS_CG(*gparams, ils_tol, ils_min_iter, ils_sigma, ils_max_iter, ils_tik_lam);
    } else if (gparams->ils_method == NLCG) {
        ils_solver = new ILS_NLCG(*gparams, ils_sigma, ils_max_iter, ils_tik_lam);
    } else if (gparams->ils_method == BiCGstabl) {
        ils_solver = new ILS_BiCGstabl(*gparams, ils_tol, ils_sigma, ils_max_iter, ils_tik_lam);
    } else {
        spdlog::error("SolverOSQP::solve()  Unknown Indirect Linear Solver method");
        return SolveResult{};
    }
    ils_solver->set_workspace(ws);

    total_Ax_size = 0;
    for (int i = 0; i < gparams->all_op.size(); i++) {
        total_Ax_size += gparams->all_op[i]->Ax_size;
    }
    r_primal.setZero(total_Ax_size);

    // ===============================================
    //  OSQP iterations
    // ===============================================
    int total_feval = 0;
    for (iiter = 0; iiter < max_iter; ++iiter) {
        spdlog::trace("Starting OSQP iteration {:d} SolverOSQP::solve", iiter);

        if (iiter > 0) {
            Xhat = ils_solver->solve(X);
        } else {
            Xhat = X;
        }

        // if (Xhat.array().isNaN().any()) {
        if ((Xhat.array().abs() > 10).any()) {
            // spdlog::error("NaN detected in Xhat at iteration {:d}. Stopping solver.", iiter);
            spdlog::error("Large values detected in Xhat at iteration {:d}. Stopping solver.", iiter);
            break;
        };

        // Update all constraints (do prox operations)
        update(Xhat);

        X = gamma_x * Xhat + (1 - gamma_x) * X;

        get_residuals(X);

        if ((logger(X) > 0) && (iiter > min_iter)) {
            break;
        }

        total_feval += ils_solver->hist_n_iter.back();
        if (total_feval > max_feval) {
            spdlog::info("Maximum function evaluations reached");
            break;
        }
    }

    SolveResult result;
    result.X = X;
    result.n_iter = iiter;
    result.dt = gparams->dt;
    final_log(X, result);

    delete ils_solver;
    spdlog::trace("Finished SolverOSQP::solve");

    return result;
}

void SolverOSQP::update(Eigen::VectorXd &X) {
    spdlog::trace("Starting SolverOSQP::update");

    double max_scale = 2.0;
    Eigen::VectorXd all_weight_scale;
    all_weight_scale.setZero(gparams->all_op.size());
    bool needs_reweight = false;

    for (int i = 0; i < gparams->all_op.size(); i++) {
        Operator *op = gparams->all_op[i].get();
        WorkspaceOSQP &w = osqp_ws[i];

        // s = Ax
        op->forward_op(X, w.s1);

        // z = prox(as + 1-a)z0 + p^-1y0)
        w.z1 = w.gamma * w.s1 + (1 - w.gamma) * w.z0 + w.y0 / w.weight;
        op->prox(w.z1);

        // y = y0 + p*(as + (1-a)z0 - z1)
        w.y1 = w.y0 + w.weight * (w.gamma * w.s1 + (1 - w.gamma) * w.z0 - w.z1);

        w.y0 = w.y1;
        w.z0 = w.z1;

        if (iiter == 2) {
            w.init_reweight();
        } else if ((iiter > 5) && (iiter % 50 == 0)) {
            w.reweight(iiter, max_scale, max_scale, 1.0e-8);
            needs_reweight = true;
        }
        all_weight_scale(i) = w.weight_scale;
    }

    if (needs_reweight) {

        std::ostringstream oss;
        oss << "ReWeights: ";
        for (int i = 0; i < gparams->all_op.size(); i++) {
            oss << all_weight_scale(i) << " ";
        }
        // spdlog::info("{:d}  {}", iiter, oss.str());

        if ((all_weight_scale.array() >= max_scale).all()) {
            // spdlog::warn("All weight scales above max_scale., scaling to max");
            // all_weight_scale = all_weight_scale / all_weight_scale.maxCoeff();
            all_weight_scale.setOnes();

            oss.str("");
            oss << "Re-ReWeights: ";
            for (int i = 0; i < gparams->all_op.size(); i++) {
                oss << all_weight_scale(i) << " ";
            }
            // spdlog::info("{:d}  {}", iiter, oss.str());
        }

        for (int i = 0; i < gparams->all_op.size(); i++) {
            Operator *op = gparams->all_op[i].get();
            WorkspaceOSQP &w = osqp_ws[i];

            if (all_weight_scale(i) > max_scale) {
                all_weight_scale(i) = max_scale;
            } else if (all_weight_scale(i) < 1.0 / max_scale) {
                all_weight_scale(i) = 1.0 / max_scale;
            }
            w.weight *= all_weight_scale(i);
        }
    }

    spdlog::trace("Finished SolverOSQP::update");
}

void SolverOSQP::get_residuals(Eigen::VectorXd &X) {

    //  Get dimensions (shouldn't be needed every iteration)
    int N_rows = 0;
    for (int i = 0; i < gparams->all_op.size(); i++) {
        N_rows += gparams->all_op[i]->Ax_size;
    }

    int N_cols = gparams->N * gparams->Naxis;

    // Calculate relevant variables for residuals
    Eigen::VectorXd Ax;
    Ax.setZero(N_rows);

    Eigen::VectorXd z;
    z.setZero(N_rows);

    Eigen::VectorXd y;
    y.setZero(N_rows);

    Eigen::VectorXd Aty;
    Aty.setZero(N_cols);

    int row_start = 0;
    for (int i = 0; i < gparams->all_op.size(); i++) {
        Operator *op = gparams->all_op[i].get();
        WorkspaceOSQP &w = osqp_ws[i];

        op->Ax_temp.setZero();
        op->forward_op(X, op->Ax_temp);
        Ax.segment(row_start, op->Ax_size) = op->Ax_temp;
        z.segment(row_start, op->Ax_size) = w.z1;
        y.segment(row_start, op->Ax_size) = w.y1;

        op->x_temp.setZero();
        op->transpose_op(w.y1, op->x_temp);
        Aty.array() += op->x_temp.array();

        row_start += op->Ax_size;
    }

    if (extra_debug) {
        debug_solver.hist_Ax.push_back(Ax);
        debug_solver.hist_z.push_back(z);
        debug_solver.hist_y.push_back(y);
        debug_solver.hist_Aty.push_back(Aty);
        debug_solver.hist_X.push_back(X);

        std::vector<double> weight_vec;
        std::vector<double> gamma_vec;
        for (int i = 0; i < gparams->all_op.size(); i++) {
            Operator *op = gparams->all_op[i].get();
            WorkspaceOSQP &w = osqp_ws[i];
            weight_vec.push_back(w.weight);
            gamma_vec.push_back(w.gamma);
        }
        debug_solver.hist_weight.push_back(weight_vec);
        debug_solver.hist_gamma.push_back(gamma_vec);

        debug_solver.hist_gamma_x.push_back(gamma_x);
    }

    // Update feasibility metrics
    for (int i = 0; i < gparams->all_op.size(); i++) {
        Operator *op = gparams->all_op[i].get();
        op->forward_op(X, op->Ax_temp);
        op->get_feas(op->Ax_temp);
        op->check(op->Ax_temp);
    }
}

} // namespace Gropt
