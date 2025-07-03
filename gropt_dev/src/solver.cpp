#include "spdlog/spdlog.h"

#include "solver.hpp"

namespace Gropt {

Solver::Solver(GroptParams &_gparams)
{
    gparams = &_gparams;
} 

void Solver::solve()
{
    spdlog::warn("Solver::solve is not implemented for the base class.");
} 

int Solver::logger(Eigen::VectorXd &X)
{
    bool do_print = (gparams->iiter % log_interval == 0);
    int all_feasible = 1;
    
    if (do_print) {
        spdlog::debug(" ");
        spdlog::debug("================= Solver Iteration {:04d} =================", gparams->iiter);
        spdlog::debug(" Last CG n_iter: {:d}   ||x|| = {:.2e}", ils_solver->hist_n_iter.back(), X.norm());
        spdlog::debug("          Name      Feasibile   Weight     Gamma     r_feas");
        spdlog::debug("------------------------------------------------------------------");
    }
    for (int i = 0; i < gparams->all_op.size(); i++) {
        if (do_print) {
            spdlog::debug("    {:^16}    {:d}       {:.1e}    {:.1e}   {:.1e}", 
                gparams->all_op[i]->name, gparams->all_op[i]->hist_feas.back(), gparams->all_op[i]->weight, gparams->all_op[i]->gamma, gparams->all_op[i]->hist_r_feas.back());
        }
        if (gparams->all_op[i]->hist_feas.back() == 0) {
            all_feasible = 0;
        }
       
    }
    
    hist_cg_iter.push_back(ils_solver->hist_n_iter.back());
    
    return all_feasible;
} 

void Solver::final_log(Eigen::VectorXd &X) 
{
    
    gparams->final_good = 1;

    spdlog::info(" ");
    spdlog::info("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ");
    spdlog::info("======================== Final Stats ========================", gparams->iiter);
    spdlog::info("  Iteration = {:d}   Total f_eval = {:d}", gparams->iiter, std::accumulate(ils_solver->hist_n_iter.begin(), ils_solver->hist_n_iter.end(), 0));
    spdlog::info("  ||x|| = {:.2e}", X.norm());
    spdlog::info(" ");
    spdlog::info("          Name      Feasibile   min(Ax)       max(Ax)      tol0 ");
    spdlog::info("-------------------------------------------------------------");
    for (int i = 0; i < gparams->all_op.size(); i++) {
        Operator *op = gparams->all_op[i];
        op->Ax_temp.setZero();
        op->forward_op(X, op->Ax_temp);
        
        spdlog::info("    {:^16}    {:d}       {: .2e}    {: .2e}    {: .2e}", 
            op->name, op->hist_feas.back(), op->Ax_temp.minCoeff()-op->target, op->Ax_temp.maxCoeff()-op->target, op->tol0);

        if (op->hist_feas.back() == 0) {
            gparams->final_good = 0;
        }
    }
}

} // namespace Gropt