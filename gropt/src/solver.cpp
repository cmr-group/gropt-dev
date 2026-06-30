#include "spdlog/spdlog.h"

#include "op_bvalue.hpp"
#include "solver.hpp"

namespace Gropt {

SolveResult Solver::solve(GroptParams &_gparams) {
    spdlog::warn("Solver::solve is not implemented for the base class.");
    return SolveResult{};
}

int Solver::logger(Eigen::VectorXd &X) {
    bool do_print = (iiter % log_interval == 0);
    int all_feasible = 1;

    if (do_print) {
        spdlog::debug(" ");
        spdlog::debug("================= Solver Iteration {:04d} =================", iiter);
        spdlog::debug(" Last CG n_iter: {:d}   ||x|| = {:.2e}", ils_solver->hist_n_iter.back(), X.norm());
        spdlog::debug("          Name      Feasibile   Weight     Gamma     r_feas");
        spdlog::debug("------------------------------------------------------------------");
    }
    for (int i = 0; i < gparams->all_op.size(); i++) {
        if (do_print) {
            spdlog::debug("    {:^16}    {:d}       {:.1e}    {:.1e}   {:.1e}", gparams->all_op[i]->name,
                          gparams->all_op[i]->hist_feas.back(), ws[i]->weight, ws[i]->gamma,
                          gparams->all_op[i]->hist_r_feas.back());
        }
        if (gparams->all_op[i]->hist_feas.back() == 0) {
            all_feasible = 0;
        }
    }

    hist_cg_iter.push_back(ils_solver->hist_n_iter.back());

    return all_feasible;
}

void Solver::final_log(Eigen::VectorXd &X, SolveResult &result) {

    result.converged = true;

    result.n_feval = std::accumulate(ils_solver->hist_n_iter.begin(), ils_solver->hist_n_iter.end(), 0);

    spdlog::info(" ");
    spdlog::info("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ");
    spdlog::info("======================== Final Stats ========================", iiter);
    spdlog::info("  Iteration = {:d}   Total f_eval = {:d}", iiter, result.n_feval);
    spdlog::info("  ||x|| = {:.2e}", X.norm());
    spdlog::info(" ");
    spdlog::info("          Name      Feasible    min(Ax)       max(Ax)      target        tol0 ");
    spdlog::info("---------------------------------------------------------------------------------");
    for (int i = 0; i < gparams->all_op.size(); i++) {
        Operator *op = gparams->all_op[i].get();
        op->Ax_temp.setZero();
        op->forward_op(X, op->Ax_temp);
        op->check(op->Ax_temp); // fresh feasibility of the RETURNED waveform (not the last loop iterate)

        // Physical (un-normalized) constraint values for the table: forward() is the raw A*X in the
        // same units as target/tol0, whereas forward_op() above is divided by spec_norm (and equil-
        // scaled) for the solver. Use a local vector so op->Ax_temp (the check() input) is untouched.
        Eigen::VectorXd Ax_phys(op->Ax_size);
        op->forward(X, Ax_phys);

        spdlog::info("    {:^16}    {:d}       {: .2e}    {: .2e}    {: .2e}    {: .2e}", op->name,
                     op->hist_feas.back(), Ax_phys.minCoeff(), Ax_phys.maxCoeff(), op->target, op->tol0);

        if (op->hist_feas.back() == 0) {
            result.converged = false;
        }
    }

    // Report the final b-value if a b-value operator exists, whether it was added as a
    // constraint (all_op) or as a maximization objective (all_obj).
    auto report_bvalue = [&](std::vector<std::unique_ptr<Operator>> &ops) {
        for (auto &op_ptr : ops) {
            Operator *op = op_ptr.get();
            if (op->name == "b-value") {
                Op_BValue *op_bvalue = dynamic_cast<Op_BValue *>(op);
                if (X.array().isNaN().any()) {
                    result.bvalue = 0;
                } else if (op_bvalue != nullptr) {
                    result.bvalue = op_bvalue->get_bvalue(X);
                } else {
                    spdlog::warn("Operator named 'b-value' is not an Op_BValue instance.");
                }
            }
        }
    };
    report_bvalue(gparams->all_op);
    report_bvalue(gparams->all_obj);
}

WarmStart Solver::capture_warmstart(const Eigen::VectorXd &X) {
    // Snapshot the current ADMM state: primal X, plus each operator's dual (ws->y1), penalty
    // (weight) and relaxation (gamma), tagged by unique_name and the operator's Ax-block layout
    // so it can be matched and resized in a later solve. The consensus z is deliberately NOT
    // stored -- it is regenerated as z = A*X when the snapshot is loaded (see warmstart.hpp).
    WarmStart w;
    if (gparams == nullptr) return w;
    w.active = true;
    w.N = gparams->N;
    w.Naxis = gparams->Naxis;
    w.dt = gparams->dt;
    w.X = X;
    w.fixer = gparams->pdata.fixer;
    for (size_t i = 0; i < gparams->all_op.size() && i < ws.size(); i++) {
        Operator *op = gparams->all_op[i].get();
        OpWarmState st;
        st.key = op->unique_name;
        st.y = ws[i]->y1;
        st.weight = ws[i]->weight;
        st.gamma = ws[i]->gamma;
        st.spec_norm = op->spec_norm; // normalization at capture, to rescale the dual on load
        st.blocks = op->Ax_block_lengths();
        w.ops.push_back(st);
    }
    return w;
}

void Solver::set_general_params(int min_iter, int max_iter, int log_interval, double gamma_x, int max_feval,
                                int obj_patience) {
    this->min_iter = min_iter;
    this->max_iter = max_iter;
    this->log_interval = log_interval;
    this->gamma_x = gamma_x;
    this->max_feval = max_feval;
    this->obj_patience = obj_patience;
}

void Solver::set_ils_params(double ils_tol, int ils_max_iter, int ils_min_iter, double ils_sigma, double ils_tik_lam) {
    this->ils_tol = ils_tol;
    this->ils_max_iter = ils_max_iter;
    this->ils_min_iter = ils_min_iter;
    this->ils_sigma = ils_sigma;
    this->ils_tik_lam = ils_tik_lam;
}

} // namespace Gropt
