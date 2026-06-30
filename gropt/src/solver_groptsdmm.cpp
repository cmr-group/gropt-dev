#include "spdlog/spdlog.h"

#include "ils.hpp"
#include "ils_bicgstabl.hpp"
#include "ils_cg.hpp"
#include "ils_nlcg.hpp"
#include "solver_groptsdmm.hpp"
#include "op_bvalue.hpp"

namespace Gropt {

SolveResult SolverGroptSDMM::solve(GroptParams &_gparams) {
    spdlog::trace("Starting SolverGroptSDMM::solve");
    gparams = &_gparams;
    if (gparams->op_prep_status != gparams->N) {
        spdlog::info("Operators do not seem prepared, calling prepare()");
        gparams->prepare();
    }

    // Resolve the starting primal.
    Eigen::VectorXd X0_init = resolve_initial_primal();
    init_workspaces(X0_init);

    // Relinearize iterate-dependent equality rows (e.g. a projected concomitant) at the ACTUAL
    // starting iterate.
    if (gparams->eq_proj_dynamic) {
        gparams->build_eq_proj(X0_init);
    }

    Eigen::VectorXd X = X0_init;
    if (gparams->eq_proj.active) {
        gparams->eq_proj.project_affine(X); // start on the equality-feasible subspace (free DOFs only)
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
        spdlog::error("SolverGroptSDMM::solve()  Unknown Indirect Linear Solver method");
        return SolveResult{};
    }
    ils_solver->set_workspace(ws);
    ils_solver->collect_debug = extra_debug; // gate the extra per-CG-iter curvature computation

    total_Ax_size = 0;
    for (int i = 0; i < gparams->all_op.size(); i++) {
        total_Ax_size += gparams->all_op[i]->Ax_size;
    }
    r_primal.setZero(total_Ax_size);

    // Locate a b-value operator (constraint or objective) for per-iteration b-value logging
    Op_BValue *bval_op = nullptr;
    if (extra_debug) {
        for (auto &op : gparams->all_op) {
            if (op->name == "b-value") bval_op = dynamic_cast<Op_BValue *>(op.get());
        }
        for (auto &op : gparams->all_obj) {
            if (op->name == "b-value") bval_op = dynamic_cast<Op_BValue *>(op.get());
        }
    }

    // The scoring objective: the optimization TARGET used for best-feasible selection (e.g. b-value
    // max). Objective-path constraint penalties (concomitant augmented Lagrangian) set is_score_obj
    // = false and are skipped here
    Operator *score_obj = nullptr;
    for (auto &o : gparams->all_obj) {
        if (o->is_score_obj) {
            score_obj = o.get();
            break;
        }
    }

    // Best feasible solution. For a single scoring objective, keep the feasible iterate that is best
    // in the direction set by the SIGN of obj_weight (<0 maximize ||A x||^2, >0 minimize).
    Eigen::VectorXd best_X;
    double best_score = 0.0; // dir * ||A_obj x||^2 of the best feasible iterate (lower = better)
    bool has_best = false;
    int iters_since_improve = 0;

    int total_feval = 0;
    for (iiter = 0; iiter < max_iter; ++iiter) {
        spdlog::trace("Starting GroptSDMM iteration {:d} SolverGroptSDMM::solve", iiter);

        if (iiter > 0) {
            // Relinearize iterate-dependent equality rows (e.g. concomitant) at the current X and
            // rebuild the projector before the solve; the CG solve then projects onto the new subspace.
            if (gparams->eq_proj_dynamic) {
                gparams->build_eq_proj(X);
            }
            // Refresh augmented-Lagrangian objective state (freeze linearization, advance dual) at X.
            for (auto &o : gparams->all_obj) {
                o->update_obj_state(X);
            }
            Xhat = ils_solver->solve(X);
        } else {
            Xhat = X;
        }

        if (((Xhat.array().abs() > 10).any()) || (Xhat.array().isNaN().any())) {
            spdlog::error("Large values detected in Xhat at iteration {:d}. Stopping solver.", iiter);
            break;
        }

        // Update all constraints (do prox operations)
        update(Xhat);

        X = gamma_x * Xhat + (1 - gamma_x) * X;

        get_residuals(X);

        if (extra_debug) {
            record_debug(X, bval_op);
        }

        int all_feasible = logger(X);
        if (extra_debug) {
            debug_solver.hist_all_feas.push_back(all_feasible);
        }
        if (all_feasible && (iiter > min_iter)) {
            if (score_obj == nullptr) {
                // No scoring objective (none, or only constraint-penalty objectives): first feasible.
                best_X = X;
                best_warmstart = capture_warmstart(X); // consistent feasible snapshot for warm starts
                has_best = true;
                debug_solver.best_feasible_iter = iiter;
                break;
            }
            // Scoring objective: score = (sign of obj_weight) * ||A_obj x||^2 (lower = better).
            Operator *obj = score_obj;
            obj->Ax_temp.setZero();
            obj->forward_op(X, obj->Ax_temp);
            double dir = (obj->obj_weight < 0.0) ? -1.0 : 1.0; // -1 maximize, +1 minimize
            double score = dir * obj->Ax_temp.squaredNorm();
            if (!has_best || score < best_score) {
                double margin = obj_rtol * (best_score < 0.0 ? -best_score : best_score);
                bool significant = !has_best || (score < best_score - margin);
                best_X = X;
                best_warmstart = capture_warmstart(X); // consistent feasible snapshot for warm starts
                best_score = score;
                has_best = true;
                debug_solver.best_feasible_iter = iiter;
                if (significant) {
                    iters_since_improve = 0;
                } else if (++iters_since_improve >= obj_patience) {
                    break;
                }
            } else if (++iters_since_improve >= obj_patience) {
                break;
            }
        }

        total_feval += ils_solver->hist_n_iter.back();
        if (total_feval > max_feval) {
            spdlog::info("Maximum function evaluations reached");
            break;
        }
    }

    // If no feasible iterate was ever found, fall back to snapshotting the final state so
    // get_warmstart() still returns something usable (flagged by the caller as infeasible-sourced).
    if (!best_warmstart.active) {
        best_warmstart = capture_warmstart(X);
    }

    SolveResult result;
    if (has_best) {
        result.X = best_X;
    } else {
        result.X = X; // no feasible iterate found; return the latest
    }
    result.n_iter = iiter;
    result.dt = gparams->dt;
    final_log(result.X, result);

    // Copy the inner-solver histories out before the solver is freed (generic to the
    // IndirectLinearSolver base, so this keeps working if the ILS method is changed).
    if (extra_debug) {
        debug_solver.hist_cg_iter = ils_solver->hist_n_iter;
        debug_solver.hist_cg_rnorm0 = ils_solver->hist_rnorm0;
        debug_solver.hist_cg_rnorm = ils_solver->hist_rnorm;
        debug_solver.hist_cg_bnorm0 = ils_solver->hist_bnorm0;
        debug_solver.hist_cg_neg_curv = ils_solver->hist_neg_curv;
        debug_solver.hist_cg_min_curv = ils_solver->hist_min_curv;
        debug_solver.hist_cg_max_curv = ils_solver->hist_max_curv;
    }

    delete ils_solver;
    spdlog::trace("Finished SolverGroptSDMM::solve");

    return result;
}

// --- solve() helpers ----------------------------------------------------------------------------- //

Eigen::VectorXd SolverGroptSDMM::resolve_initial_primal() {
    // With a warm start, resize the snapshot's waveform onto this problem's grid: free segments
    // interpolate, fixed segments come from set_vals (see warmstart.hpp). Fall back to the cold X0 if
    // the snapshot is incompatible (different Naxis or a different free/fixed segment topology -- a
    // structural change, not a resize).
    Eigen::VectorXd X0_init = gparams->pdata.X0;
    if (warmstart.active) {
        bool compatible = (warmstart.Naxis == gparams->Naxis) &&
                          (ws_free_run_count(warmstart.fixer, warmstart.Naxis) ==
                           ws_free_run_count(gparams->pdata.fixer, gparams->Naxis));
        if (compatible) {
            X0_init = ws_resize_waveform(warmstart.X, warmstart.fixer, gparams->pdata.fixer,
                                         gparams->pdata.set_vals, gparams->Naxis);
        } else {
            spdlog::warn("Warm start incompatible (Naxis or free-segment count mismatch); "
                         "ignoring it and starting cold.");
            warmstart.active = false;
        }
    }
    return X0_init;
}

void SolverGroptSDMM::init_workspaces(const Eigen::VectorXd &X0_init) {
    sdmm_ws.resize(gparams->all_op.size());
    for (int i = 0; i < gparams->all_op.size(); i++) {
        Operator *op = gparams->all_op[i].get();

        // Initial ADMM weight, scaled by the operator's problem-specific weight_mod.
        sdmm_ws[i].weight = 1.0 * op->weight_mod;

        // Pin the relaxation gamma for operators that ask for it (e.g. nonconvex prox), so the
        // reweighter's over-relaxation can't destabilize them. The weight still adapts independently.
        if (op->fix_gamma) {
            sdmm_ws[i].gamma = op->gamma_fix;
            sdmm_ws[i].do_gamma = false;
        }

        sdmm_ws[i].init(op->Ax_size);
        sdmm_ws[i].prep(*op, X0_init); // z = A*X0_init  (regenerated from the primal, cold or warm)

        // Warm-start dual injection: match by unique_name, resize the dual y onto this operator's
        // grid, and seed its penalty/relaxation. An unmatched operator (e.g. a newly added
        // constraint) stays cold (y = 0) so the early iterations discover it -- a constraint-space
        // homotopy. The reweighter re-adapts the seeded weight from here.
        if (warmstart.active) {
            const OpWarmState *st = warmstart.find(op->unique_name);
            if (st != nullptr) {
                sdmm_ws[i].y0 = ws_resize_dual(*st, op->Ax_block_lengths(), op->spec_norm);
                sdmm_ws[i].y1 = sdmm_ws[i].y0;
                sdmm_ws[i].weight = st->weight;
                if (!op->fix_gamma) sdmm_ws[i].gamma = st->gamma;
            }
        }
    }

    // One-shot: the loaded warm start applies to exactly this solve. Reset it so a later solve() on
    // the same solver is cold.
    warmstart = WarmStart{};

    // Populate the base-class ws pointers into the typed SDMM workspaces.
    ws.resize(sdmm_ws.size());
    for (int i = 0; i < sdmm_ws.size(); i++) {
        ws[i] = &sdmm_ws[i];
    }
}

void SolverGroptSDMM::record_debug(Eigen::VectorXd &X, Op_BValue *bval_op) {
    Eigen::VectorXd Ax(total_Ax_size), z(total_Ax_size), y(total_Ax_size);
    Eigen::VectorXd Aty = Eigen::VectorXd::Zero(X.size());
    std::vector<double> weight_vec, gamma_vec, rfeas_vec;
    std::vector<double> con_pull_vec; // per-op ||Aᵀy||: this op's pull on the x-subproblem RHS
    std::vector<int> feas_vec;
    int row_start = 0;
    for (int i = 0; i < gparams->all_op.size(); i++) {
        Operator *op = gparams->all_op[i].get();
        WorkspaceSDMM &w = sdmm_ws[i];

        op->Ax_temp.setZero();
        op->forward_op(X, op->Ax_temp);
        Ax.segment(row_start, op->Ax_size) = op->Ax_temp;
        z.segment(row_start, op->Ax_size) = w.z1;
        y.segment(row_start, op->Ax_size) = w.y1;

        op->x_temp.setZero();
        op->transpose_op(w.y1, op->x_temp);
        Aty.array() += op->x_temp.array();
        con_pull_vec.push_back(op->x_temp.norm());

        weight_vec.push_back(w.weight);
        gamma_vec.push_back(w.gamma);
        rfeas_vec.push_back(op->hist_r_feas.back());
        feas_vec.push_back(op->hist_feas.back());

        row_start += op->Ax_size;
    }

    // --- Balance diagnostics: objective (DCA/RHS) pull vs constraint pulls ---
    // ||g_obj|| = ||Σ_obj -obj_weight·AᵀA·X||, the frozen pull add_obj_rhs injects into the
    // x-subproblem RHS each solve. Compared against ||Σ Aᵀy|| (total constraint pull). Both are in
    // normalized x-space (forward_op/transpose_op already divide by spec_norm), so the ratio is
    // apples-to-apples.
    Eigen::VectorXd obj_pull = Eigen::VectorXd::Zero(X.size());
    for (auto &obj : gparams->all_obj) {
        obj->Ax_temp.setZero();
        obj->x_temp.setZero();
        obj->forward_op(X, obj->Ax_temp);
        obj->transpose_op(obj->Ax_temp, obj->x_temp);
        obj_pull.array() += -obj->obj_weight * obj->x_temp.array();
    }
    debug_solver.hist_obj_pull.push_back(obj_pull.norm());
    debug_solver.hist_con_pull.push_back(Aty.norm());
    debug_solver.hist_con_pull_op.push_back(con_pull_vec);

    debug_solver.hist_X.push_back(X);
    debug_solver.hist_Ax.push_back(Ax);
    debug_solver.hist_z.push_back(z);
    debug_solver.hist_y.push_back(y);
    debug_solver.hist_Aty.push_back(Aty);
    debug_solver.hist_weight.push_back(weight_vec);
    debug_solver.hist_gamma.push_back(gamma_vec);
    debug_solver.hist_gamma_x.push_back(gamma_x);
    debug_solver.hist_r_feas.push_back(rfeas_vec);
    debug_solver.hist_feas.push_back(feas_vec);
    if (bval_op != nullptr) {
        debug_solver.hist_bvalue.push_back(bval_op->get_bvalue(X));
    }
}

void SolverGroptSDMM::update(Eigen::VectorXd &X) {
    spdlog::trace("Starting SolverGroptSDMM::update");

    // Per-iteration ADMM residual traces (only populated when extra_debug is on)
    std::vector<double> dbg_r_prim;
    std::vector<double> dbg_r_dual;

    for (int i = 0; i < gparams->all_op.size(); i++) {
        Operator *op = gparams->all_op[i].get();
        if (op->use_projection) continue; // moments handled by null-space projection, not ADMM
        WorkspaceSDMM &w = sdmm_ws[i];

        // s = Ax
        op->forward_op(X, w.s1);

        // z = prox(as + 1-a)z0 + p^-1y0)
        w.z1 = w.gamma * w.s1 + (1 - w.gamma) * w.z0 + w.y0 / w.weight;
        op->admm_weight = w.weight; // expose rho to penalty proxes (Op_TV divides by it); box proxes ignore
        op->prox(w.z1);

        // y = y0 + p*(as + (1-a)z0 - z1)
        w.y1 = w.y0 + w.weight * (w.gamma * w.s1 + (1 - w.gamma) * w.z0 - w.z1);

        // Boyd ADMM diagnostic residuals.
        if (extra_debug) {
            Eigen::VectorXd dz = w.z1 - w.z0;
            op->transpose_op(dz, op->x_temp); // x_temp is x-space scratch, recomputed in get_residuals
            dbg_r_prim.push_back((w.s1 - w.z1).norm());
            dbg_r_dual.push_back(w.weight * op->x_temp.norm());
        }

        if ((w.do_rw) && (iiter > rw_interval) && (iiter % rw_interval == 0)) {
            w.reweight(rw_eps, rw_e_corr, rw_scalelim);
        }

        w.y0 = w.y1;
        w.z0 = w.z1;
    }

    if (extra_debug) {
        debug_solver.hist_r_prim.push_back(dbg_r_prim);
        debug_solver.hist_r_dual.push_back(dbg_r_dual);
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
