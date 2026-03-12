#include "spdlog/spdlog.h"

#include "ils.hpp"
#include "ils_bicgstabl.hpp"
#include "ils_cg.hpp"
#include "ils_nlcg.hpp"
#include "solver_groptsdmm.hpp"

namespace Gropt {

SolveResult SolverGroptSDMM::solve(GroptParams &_gparams) {
    spdlog::trace("Starting SolverGroptSDMM::solve");
    gparams = &_gparams;
    if (gparams->op_prep_status != gparams->N) {
        spdlog::info("Operators do not seem prepared, calling prepare()");
        gparams->prepare();
    }

    // Initialize per-operator SDMM workspaces
    sdmm_ws.resize(gparams->all_op.size());
    for (int i = 0; i < gparams->all_op.size(); i++) {
        Operator *op = gparams->all_op[i].get();

        // Set initial weight based on operator type
        sdmm_ws[i].weight = 1.0;
        // Slew, moment, bvalue, SAFE, TV operators start with higher weight
        if (op->name == "Slew" || op->name == "Moment" || op->name == "b-value" || op->name == "SAFE" ||
            op->name == "TotalVariation") {
            sdmm_ws[i].weight = 1e4;
        }
        sdmm_ws[i].weight *= op->weight_mod;

        sdmm_ws[i].init(op->Ax_size);
        sdmm_ws[i].prep(*op, gparams->pdata.X0);
    }

    // Populate base class ws pointers
    ws.resize(sdmm_ws.size());
    for (int i = 0; i < sdmm_ws.size(); i++) {
        ws[i] = &sdmm_ws[i];
    }

    Eigen::VectorXd X = gparams->pdata.X0;
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
        spdlog::error("SolverGroptSDMM::solve()  Unknown Indirect Linear Solver method");
        return SolveResult{};
    }
    ils_solver->set_workspace(ws);

    total_Ax_size = 0;
    for (int i = 0; i < gparams->all_op.size(); i++) {
        total_Ax_size += gparams->all_op[i]->Ax_size;
    }
    r_primal.setZero(total_Ax_size);

    int total_feval = 0;
    for (iiter = 0; iiter < max_iter; ++iiter) {
        spdlog::trace("Starting GroptSDMM iteration {:d} SolverGroptSDMM::solve", iiter);

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

        if (extra_debug) {
            debug_solver.hist_X.push_back(X);
        }

        if ((logger(X) > 0) && (iiter > min_iter)) {
            if (extra_iters > 0) {
                spdlog::info("First dsolved at iiter {:d}, now {:d} extra iterations", iiter, extra_iters);
                min_iter = iiter + extra_iters;
                extra_iters = 0;
            } else {
                break;
            }
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
    spdlog::trace("Finished SolverGroptSDMM::solve");

    return result;
}

void SolverGroptSDMM::update(Eigen::VectorXd &X) {
    spdlog::trace("Starting SolverGroptSDMM::update");

    for (int i = 0; i < gparams->all_op.size(); i++) {
        Operator *op = gparams->all_op[i].get();
        WorkspaceSDMM &w = sdmm_ws[i];

        // s = Ax
        op->forward_op(X, w.s1);

        // z = prox(as + 1-a)z0 + p^-1y0)
        w.z1 = w.gamma * w.s1 + (1 - w.gamma) * w.z0 + w.y0 / w.weight;
        op->prox(w.z1);

        // y = y0 + p*(as + (1-a)z0 - z1)
        w.y1 = w.y0 + w.weight * (w.gamma * w.s1 + (1 - w.gamma) * w.z0 - w.z1);

        if ((w.do_rw) && (iiter > rw_interval) && (iiter % rw_interval == 0)) {
            w.reweight(rw_eps, rw_e_corr, rw_scalelim);
        }

        w.y0 = w.y1;
        w.z0 = w.z1;
    }

    spdlog::trace("Finished SolverGroptSDMM::update");
}

void SolverGroptSDMM::get_residuals(Eigen::VectorXd &X) {
    // Update feasibility metrics
    for (int i = 0; i < gparams->all_op.size(); i++) {
        Operator *op = gparams->all_op[i].get();
        op->forward_op(X, op->Ax_temp);
        op->get_feas(op->Ax_temp);
        op->check(op->Ax_temp);
    }

    if (iiter > 2 * grw_min_infeasible && iiter % grw_interval == 0) {
        double max_feas = 0.0;
        int max_index = -1;
        for (int i = 0; i < gparams->all_op.size(); i++) {
            if (std::accumulate(gparams->all_op[i]->hist_feas.end() - grw_min_infeasible,
                                gparams->all_op[i]->hist_feas.end(), 0) == 0) {
                if (gparams->all_op[i]->hist_r_feas.back() > max_feas) {
                    max_feas = gparams->all_op[i]->hist_r_feas.back();
                    max_index = i;
                }
            }
        }
        if (max_index >= 0) {
            sdmm_ws[max_index].weight *= grw_mod;
        }
    }
}

void SolverGroptSDMM::set_sdmm_params(int rw_interval, double rw_e_corr, double rw_eps, double rw_scalelim,
                                      int grw_min_infeasible, int grw_interval, double grw_mod) {
    this->rw_interval = rw_interval;
    this->rw_e_corr = rw_e_corr;
    this->rw_eps = rw_eps;
    this->rw_scalelim = rw_scalelim;

    this->grw_min_infeasible = grw_min_infeasible;
    this->grw_interval = grw_interval;
    this->grw_mod = grw_mod;
}

} // namespace Gropt
