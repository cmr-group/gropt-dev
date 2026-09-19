#include <cmath>
#include <stdexcept>

#include "spdlog/spdlog.h"

#include "ils.hpp"
#include "ils_bicgstabl.hpp"
#include "ils_cg.hpp"
#include "ils_nlcg.hpp"
#include "solver_groptsdmm.hpp"
#include "op_bvalue.hpp"
#include "fft_tools.hpp"

namespace Gropt {

SolveResult SolverGroptSDMM::solve(GroptParams &_gparams) {
    spdlog::trace("Starting SolverGroptSDMM::solve");
    gparams = &_gparams;
    if (gparams->needs_prepare()) {
        spdlog::info("Operators do not seem prepared, calling prepare()");
        gparams->prepare();
    }

    // Only ILS_CG applies the equality projection; other inner solvers see projected constraints only
    // through the outer re-projection. Checked before the warm start is consumed.
    if (gparams->eq_proj.active && gparams->ils_method != CG) {
        const char *ils_name = (gparams->ils_method == NLCG) ? "NLCG" : "BiCGstabl";
        if (!reproject_iterate) {
            throw std::invalid_argument(std::string(ils_name) +
                                        " does not apply the equality projection; with reproject_iterate=false "
                                        "the projected constraints would be ignored. Use CG.");
        }
        spdlog::warn("{} does not apply the equality projection; projected constraints are only enforced by "
                     "reproject_iterate (slower). Use CG.", ils_name);
    }

    // Per-solve outputs and gates start fresh, so a reused solver or gparams carries nothing over.
    best_warmstart = WarmStart{};
    debug_solver = DebugSolver{};
    for (auto &o : gparams->all_obj) {
        o->obj_gate = 1.0;
    }

    Eigen::VectorXd X0_init = resolve_initial_primal();
    init_workspaces(X0_init);

    // Relinearize iterate-dependent equality rows (e.g. projected concomitant) at the starting iterate.
    if (gparams->eq_proj_dynamic) {
        gparams->build_eq_proj(X0_init);
    }

    Eigen::VectorXd X = X0_init;
    if (gparams->eq_proj.active) {
        gparams->eq_proj.project_affine(X); // start on the equality-feasible subspace (free DOFs only)
    }
    Eigen::VectorXd Xhat;

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

    LowFreqProjector lowfreq; // inactive unless cutoff_freq > 0
    if (cutoff_freq > 0.0) {
        lowfreq.setup(gparams->N, gparams->Naxis, gparams->dt, cutoff_freq, gparams->pdata.fixer,
                      cutoff_trans);
    }

    // b-value operator (constraint or objective) for the per-iteration debug history
    Op_BValue *bval_op = nullptr;
    if (extra_debug) {
        for (auto &op : gparams->all_op) {
            if (op->name == "b-value") bval_op = dynamic_cast<Op_BValue *>(op.get());
        }
        for (auto &op : gparams->all_obj) {
            if (op->name == "b-value") bval_op = dynamic_cast<Op_BValue *>(op.get());
        }
    }

    // The first objective (e.g. b-value max) scores feasible iterates for best-feasible selection.
    Operator *score_obj = gparams->all_obj.empty() ? nullptr : gparams->all_obj.front().get();

    // Best feasible iterate; obj_weight < 0 maximizes ||A x||^2, > 0 minimizes it.
    Eigen::VectorXd best_X;
    double best_score = 0.0; // dir * ||A_obj x||^2 of the best feasible iterate (lower = better)
    bool has_best = false;
    int iters_since_improve = 0;

    // Trust region: the proximal sigma grows on rejected steps and decays back toward ils_sigma.
    std::unique_ptr<StepMonitor> step_monitor;
    if (tr_enable) {
        step_monitor = make_step_monitor(tr_monitor, tr_tol, tr_bump);
        if (!step_monitor && tr_monitor != "none") {
            spdlog::warn("Unknown tr_monitor '{}'; trust region disabled.", tr_monitor);
        }
    }
    double tr_sigma = ils_sigma;

    int total_feval = 0;
    for (iiter = 0; iiter < max_iter; ++iiter) {
        spdlog::trace("Starting GroptSDMM iteration {:d} SolverGroptSDMM::solve", iiter);

        if (iiter > 0) {
            if (gparams->eq_proj_dynamic) {
                gparams->build_eq_proj(X);
            }
            if (obj_gate_enable) {
                double viol = 0.0;
                for (auto &op : gparams->all_op) {
                    viol += op->constraint_violation(X);
                }
                double g = std::exp(-viol / std::max(obj_gate_scale, 1e-30));
                for (auto &o : gparams->all_obj) {
                    o->obj_gate = g;
                }
            }

            // Inner solve with nonlinear operators linearized (frozen) at X.
            if (step_monitor) {
                // Re-solve from X with a larger sigma until the monitor accepts (or tr_max_reject).
                int nrej = 0;
                while (true) {
                    for (auto &op : gparams->all_op) op->freeze_linearization(X);
                    ils_solver->sigma = tr_sigma;
                    Xhat = ils_solver->solve(X);
                    total_feval += ils_solver->hist_n_iter.back(); // every re-solve counts toward max_feval
                    StepDecision d = step_monitor->check(*gparams, X, Xhat, tr_sigma);
                    for (auto &op : gparams->all_op) op->unfreeze_linearization();
                    if (d.accept) {
                        tr_sigma = std::max(ils_sigma, tr_sigma * tr_decay);
                        break;
                    }
                    spdlog::debug("iter {:d}: {} rejected step (signal {:.3e}, tol {:.3e}), sigma {:.2e}", iiter,
                                  step_monitor->name(), d.signal, step_monitor->tol, tr_sigma);
                    if (++nrej > tr_max_reject) {
                        spdlog::debug("iter {:d}: {} rejects exceeded, taking the damped step", iiter, nrej);
                        break; // give up: take the (most-damped) step rather than stall
                    }
                    tr_sigma *= d.sigma_scale;
                }
            } else {
                for (auto &op : gparams->all_op) op->freeze_linearization(X);
                Xhat = ils_solver->solve(X);
                total_feval += ils_solver->hist_n_iter.back();
                for (auto &op : gparams->all_op) op->unfreeze_linearization();
            }
        } else {
            Xhat = X;
        }

        if (((Xhat.array().abs() > 10).any()) || (Xhat.array().isNaN().any())) {
            spdlog::error("Large values detected in Xhat at iteration {:d}. Stopping solver.", iiter);
            break;
        }

        update(Xhat); // ADMM z/y updates (prox)

        X = gamma_x * Xhat + (1 - gamma_x) * X;

        if (reproject_iterate && gparams->eq_proj.active) {
            gparams->eq_proj.project_affine(X);
        }

        // Low-pass the iterate until cutoff_iter (< 0 = always); re-project, since it changes the moments.
        if (lowfreq.active() && (cutoff_iter < 0 || iiter < cutoff_iter)) {
            lowfreq.project(X);
            if (gparams->eq_proj.active) {
                gparams->eq_proj.project_affine(X);
            }
        }

        get_residuals(X);

        if (extra_debug) {
            record_debug(X, bval_op);
        }

        int all_feasible = logger(X);
        if (extra_debug) {
            debug_solver.hist_all_feas.push_back(all_feasible);
        }
        if (all_feasible && (iiter > min_iter)) {
            auto save_best = [&]() {
                best_X = X;
                best_warmstart = capture_warmstart(X); // consistent feasible snapshot for warm starts
                has_best = true;
                debug_solver.best_feasible_iter = iiter;
            };
            if (score_obj == nullptr) {
                save_best(); // no objective: return the first feasible iterate
                break;
            }
            // score = (sign of obj_weight) * ||A_obj x||^2, lower = better
            score_obj->Ax_temp.setZero();
            score_obj->forward_op(X, score_obj->Ax_temp);
            const double dir = (score_obj->obj_weight < 0.0) ? -1.0 : 1.0; // -1 maximize, +1 minimize
            const double score = dir * score_obj->Ax_temp.squaredNorm();
            const bool significant = !has_best || (score < best_score - obj_rtol * std::abs(best_score));
            if (!has_best || score < best_score) {
                save_best();
                best_score = score;
            }
            if (significant) {
                iters_since_improve = 0;
            } else if (++iters_since_improve >= obj_patience) {
                break;
            }
        }

        if (total_feval > max_feval) {
            spdlog::info("Maximum function evaluations reached");
            break;
        }
    }

    // No feasible iterate: snapshot the final state so get_warmstart() still returns something usable.
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

    // Copy the inner-solver histories before ils_solver is freed.
    if (extra_debug) {
        debug_solver.hist_cg_iter = ils_solver->hist_n_iter;
        debug_solver.hist_cg_rnorm0 = ils_solver->hist_rnorm0;
        debug_solver.hist_cg_rnorm = ils_solver->hist_rnorm;
        debug_solver.hist_cg_bnorm0 = ils_solver->hist_bnorm0;
    }

    delete ils_solver;
    spdlog::trace("Finished SolverGroptSDMM::solve");

    return result;
}

// --- solve() helpers ----------------------------------------------------------------------------- //

Eigen::VectorXd SolverGroptSDMM::resolve_initial_primal() {
    // Warm start: resize the snapshot waveform onto this grid (see warmstart.hpp). Start cold if its
    // axes or free/fixed layout don't match.
    Eigen::VectorXd X0_init = gparams->pdata.X0;
    if (warmstart.active) {
        bool compatible = (warmstart.Naxis == gparams->Naxis) && (warmstart.Naxis > 0) &&
                          (warmstart.X.size() == warmstart.fixer.size()) &&
                          (warmstart.fixer.size() % warmstart.Naxis == 0) &&
                          (ws_free_run_counts(warmstart.fixer, warmstart.Naxis) ==
                           ws_free_run_counts(gparams->pdata.fixer, gparams->Naxis));
        if (compatible) {
            X0_init = ws_resize_waveform(warmstart.X, warmstart.fixer, gparams->pdata.fixer,
                                         gparams->pdata.set_vals, gparams->Naxis);
        } else {
            spdlog::warn("Warm start incompatible (Naxis, size, or free-segment count mismatch); "
                         "ignoring it and starting cold.");
            warmstart.active = false;
        }
    }
    return X0_init;
}

void SolverGroptSDMM::init_workspaces(const Eigen::VectorXd &X0_init) {
    sdmm_ws.assign(gparams->all_op.size(), WorkspaceSDMM{}); // fresh: no gamma/flags from a prior solve
    for (int i = 0; i < gparams->all_op.size(); i++) {
        Operator *op = gparams->all_op[i].get();

        sdmm_ws[i].weight = 1.0 * op->weight_mod;

        sdmm_ws[i].init(op->Ax_size);
        sdmm_ws[i].prep(*op, X0_init); // z = A*X0_init (also for warm starts)

        // Warm start: seed y, weight and gamma from the snapshot operator with the same unique_name;
        // unmatched operators start cold.
        if (warmstart.active) {
            const OpWarmState *st = warmstart.find(op->unique_name);
            if (st != nullptr) {
                Eigen::VectorXd y_warm = ws_resize_dual(*st, op->Ax_block_lengths(), op->spec_norm);
                if (y_warm.size() == op->Ax_size) {
                    sdmm_ws[i].y0 = y_warm;
                    sdmm_ws[i].y1 = y_warm;
                    sdmm_ws[i].weight = st->weight;
                    sdmm_ws[i].gamma = st->gamma;
                } else {
                    spdlog::warn("Warm start for '{}' has an incompatible dual layout; starting it cold.",
                                 op->unique_name);
                }
            }
        }
    }

    // A loaded warm start applies to this solve only.
    warmstart = WarmStart{};

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

    // Linearized-objective RHS pull at X (with normalize_obj and the current obj_gate) vs ||Σ Aᵀy||.
    Eigen::VectorXd obj_pull = Eigen::VectorXd::Zero(X.size());
    for (auto &obj : gparams->all_obj) {
        obj->add_obj_rhs(X, obj_pull, gparams->normalize_obj);
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
        if (op->use_projection) continue; // enforced by the equality projection, not ADMM
        WorkspaceSDMM &w = sdmm_ws[i];

        // s = A x
        op->forward_op(X, w.s1);

        // z = prox(gamma*s + (1-gamma)*z0 + y0/rho)
        w.z1 = w.gamma * w.s1 + (1 - w.gamma) * w.z0 + w.y0 / w.weight;
        op->admm_weight = w.weight; // penalty proxes (Op_TV) need rho; constraint proxes ignore it
        op->prox(w.z1);

        // y = y0 + rho*(gamma*s + (1-gamma)*z0 - z1)
        w.y1 = w.y0 + w.weight * (w.gamma * w.s1 + (1 - w.gamma) * w.z0 - w.z1);

        // Boyd ADMM diagnostic residuals.
        if (extra_debug) {
            Eigen::VectorXd dz = w.z1 - w.z0;
            op->transpose_op(dz, op->x_temp); // x_temp is x-space scratch, recomputed in get_residuals
            dbg_r_prim.push_back((w.s1 - w.z1).norm());
            dbg_r_dual.push_back(w.weight * op->x_temp.norm());
        }

        if (bb_enable && w.do_rw && (iiter > rw_interval) && (iiter % rw_interval == 0)) {
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

    if (grw_enable && iiter > 2 * grw_min_infeasible && iiter % grw_interval == 0) {
        double max_feas = 0.0;
        int max_index = -1;
        for (int i = 0; i < gparams->all_op.size(); i++) {
            if (gparams->all_op[i]->use_projection) continue; // ADMM weight unused
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
            if (grw_balanced && grw_mod > 0.0) {
                // Keep the geometric mean of the ADMM (non-projected) constraint weights fixed.
                int K = 0;
                for (auto &op : gparams->all_op) {
                    if (!op->use_projection) K++;
                }
                if (K > 0) {
                    double renorm = std::pow(grw_mod, 1.0 / K);
                    for (int i = 0; i < static_cast<int>(gparams->all_op.size()); i++) {
                        if (!gparams->all_op[i]->use_projection) {
                            sdmm_ws[i].weight /= renorm;
                        }
                    }
                }
            }
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
