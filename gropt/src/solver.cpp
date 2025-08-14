#include "spdlog/spdlog.h"

#include "solver.hpp"

namespace Gropt {

void Solver::solve(GroptParams &_gparams)
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

    gparams->final_n_feval = std::accumulate(ils_solver->hist_n_iter.begin(), ils_solver->hist_n_iter.end(), 0);

    spdlog::info(" ");
    spdlog::info("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ");
    spdlog::info("======================== Final Stats ========================", gparams->iiter);
    spdlog::info("  Iteration = {:d}   Total f_eval = {:d}", gparams->iiter, gparams->final_n_feval);
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

void Solver::set_general_params(int min_iter, int max_iter, int log_interval, double gamma_x, int max_feval)
{
    this->min_iter = min_iter;
    this->max_iter = max_iter;
    this->log_interval = log_interval;
    this->gamma_x = gamma_x;
    this->max_feval = max_feval;
}

void Solver::set_ils_params(double ils_tol, int ils_max_iter, int ils_min_iter, double ils_sigma, double ils_tik_lam)
{
    this->ils_tol = ils_tol;
    this->ils_max_iter = ils_max_iter;
    this->ils_min_iter = ils_min_iter;
    this->ils_sigma = ils_sigma;
    this->ils_tik_lam = ils_tik_lam;
}

} // namespace Gropt