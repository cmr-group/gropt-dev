#include "spdlog/spdlog.h"

#include "ils.hpp"
#include "ils_cg.hpp"
#include "ils_nlcg.hpp"
#include "ils_bicgstabl.hpp"
#include "solver_groptsdmm.hpp"

namespace Gropt {

void SolverGroptSDMM::solve()  
{
    spdlog::trace("Starting SolverGroptSDMM::solve");

    Eigen::VectorXd X = gparams->X0;
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
        return;
    }

    total_Ax_size = 0;
    for (int i = 0; i < gparams->all_op.size(); i++) {
        total_Ax_size += gparams->all_op[i]->Ax_size;
    }
    r_primal.setZero(total_Ax_size);
    
    int iiter;
    for (iiter = 0; iiter < max_iter; ++iiter) {
        spdlog::trace("Starting GroptSDMM iteration {:d} SolverGroptSDMM::solve", iiter);
        gparams->iiter = iiter;

        if (iiter > 0) {
            Xhat = ils_solver->solve(X);
        } else {
            Xhat = X;
        }

        if (Xhat.array().isNaN().any()) {
            spdlog::error("NaN detected in Xhat at iteration {:d}. Stopping solver.", iiter);
            break;
        };
        // Xhat.array() *= gparams->fixer.array();
 
        // Update all constraints (do prox operations)        
        update(Xhat);

        X = gamma_x * Xhat + (1 - gamma_x) * X;

        get_residuals(X);

        if (extra_debug) {hist_X.push_back(X);}
        
        if ((logger(X) > 0) && (iiter > min_iter)) {break;}
    }

    gparams->final_X = X;
    final_log(X);

    delete ils_solver; 
    spdlog::trace("Finished Starting SolverGroptSDMM::solve");

}

void SolverGroptSDMM::update(Eigen::VectorXd &X)
{
    spdlog::trace("Starting SolverGroptSDMM::update");

    for (int i = 0; i < gparams->all_op.size(); i++) {
        Operator *op = gparams->all_op[i];
        
        // s = Ax
        op->forward_op(X, op->s1);

        // z = prox(as + 1-a)z0 + p^-1y0)
        op->z1 = op->gamma * op->s1 + (1 - op->gamma) * op->z0 + op->y0 / op->weight;
        op->prox(op->z1);

        // y = y0 + p*(as + (1-a)z0 - z1)
        op->y1 = op->y0 + op->weight * (op->gamma * op->s1 + (1 - op->gamma) * op->z0 - op->z1);

        if ((op->do_rw) && (gparams->iiter > op->rw_interval) && (gparams->iiter%op->rw_interval == 0)) {
            op->reweight_parsdmm();
        }

        // op->dyk = op->yk1 - op->yk0;
        // op->dzk = op->zk1 - op->zk0;

        op->y0 = op->y1;
        op->z0 = op->z1;
    }       

    spdlog::trace("Finished SolverGroptSDMM::update");
}

void SolverGroptSDMM::get_residuals(Eigen::VectorXd &X)
{
    // Update feasibility metrics
    for (int i = 0; i < gparams->all_op.size(); i++) {
        Operator *op = gparams->all_op[i];
        op->forward_op(X, op->Ax_temp);
        op->get_feas(op->Ax_temp);
        op->check(op->Ax_temp);
    }

    int N_worst_update = 20;
    if (gparams->iiter > 2*N_worst_update && gparams->iiter % N_worst_update == 0) {
        double max_feas = 0.0;
        int max_index = -1;
        for (int i = 0; i < gparams->all_op.size(); i++) {
            if (std::accumulate(gparams->all_op[i]->hist_feas.end()-N_worst_update, gparams->all_op[i]->hist_feas.end(), 0) == 0) {
                if (gparams->all_op[i]->hist_r_feas.back() > max_feas) {
                    max_feas = gparams->all_op[i]->hist_r_feas.back();
                    max_index = i;
                }
            } 
        }
        if (max_index >= 0) {
            gparams->all_op[max_index]->weight *= 1.5;
        }
    }


}

/*
void SolverGroptSDMM::get_residuals(Eigen::VectorXd &X)
{
    double max_Ax = 0.0;
    double max_z = 0.0;
    double max_Aty = 0.0;
    double max_Px = 0.0;
    double max_q = 0.0;
    
    Px.setZero();
    r_dual.setZero();
    for (int i = 0; i < gparams->all_obj.size(); i++) {
        Operator *op = gparams->all_obj[i];
        op->x_temp_obj.setZero();
        op->add_obj(X, op->x_temp_obj);
        Px += op->x_temp_obj;
    }
    max_Px = Px.cwiseAbs().maxCoeff();
    r_dual += Px;


    int i_total = 0;
    for (int i = 0; i < gparams->all_op.size(); i++) {
        Operator *op = gparams->all_op[i];
        
        // ---  r_primal = Ax - z
        op->forward_op(X, op->Ax_temp);
        // op->r_primal = op->Ax_temp - op->zk1;
        r_primal.segment(i_total, op->Ax_size) = op->Ax_temp - op->z1;
        
        max_Ax = std::max(op->Ax_temp.cwiseAbs().maxCoeff(), max_Ax);
        max_z = std::max(op->z1.cwiseAbs().maxCoeff(), max_z);


        // ---  r_dual = Px + q + A'y
        op->transpose_op(op->y1, op->x_temp);
        // op->r_dual = Px + op->x_temp; 
        r_dual += op->x_temp;

        max_Aty = std::max(op->x_temp.cwiseAbs().maxCoeff(), max_Aty);
        
        i_total += op->Ax_size;
    }       

    double max_primal = std::max({max_Ax, max_z});
    double max_dual = std::max({max_Px, max_q, max_Aty});

    double scale_numer = r_primal.cwiseAbs().maxCoeff() / max_primal;
    double scale_denom = r_dual.cwiseAbs().maxCoeff() / max_dual;
    double scale_factor = sqrt(scale_numer / scale_denom);

    // if (gparams->iiter > 0 && gparams->iiter % 10 == 0 && (scale_factor > 5.0 || scale_factor < 0.2)) {
    //     for (int i = 0; i < gparams->all_op.size(); i++) {
    //         gparams->all_op[i]->weight *= scale_factor;
    //     }
    // }

    // spdlog::debug("~~ max_Px {:.2e}  scale_factor {:.2e}", max_Px, scale_factor);

}
*/

} // namespace Gropt