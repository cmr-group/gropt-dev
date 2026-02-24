#include "spdlog/spdlog.h"

#include "ils.hpp"
#include "ils_cg.hpp"
#include "ils_nlcg.hpp"
#include "ils_bicgstabl.hpp"
#include "solver_osqp.hpp"

namespace Gropt
{

    SolveResult SolverOSQP::solve(GroptParams &_gparams)
    {
        spdlog::trace("Starting SolverOSQP::solve");
        gparams = &_gparams;
        if (gparams->op_prep_status != gparams->N)
        {
            spdlog::info("Operators do not seem prepared, calling prepare()");
            gparams->prepare();
        }

        // Initialize per-operator OSQP workspaces
        osqp_ws.resize(gparams->all_op.size());
        for (int i = 0; i < gparams->all_op.size(); i++)
        {
            Operator *op = gparams->all_op[i].get();

            // Set initial weight based on operator type
            osqp_ws[i].weight = 1.0;
            // Slew, moment, bvalue, SAFE, TV operators start with higher weight
            if (op->name == "Slew" || op->name == "Moment" || op->name == "b-value" ||
                op->name == "SAFE" || op->name == "TotalVariation")
            {
                osqp_ws[i].weight = 1e4;
            }
            osqp_ws[i].weight *= op->weight_mod;

            osqp_ws[i].init(op->Ax_size);
            osqp_ws[i].prep(*op, gparams->pdata.X0);
        }

        // Populate base class ws pointers
        ws.resize(osqp_ws.size());
        for (int i = 0; i < osqp_ws.size(); i++)
        {
            ws[i] = &osqp_ws[i];
        }

        Eigen::VectorXd X = gparams->pdata.X0;
        Eigen::VectorXd Xhat;

        Px.setZero(X.size());
        r_dual.setZero(X.size());

        if (gparams->ils_method == CG)
        {
            ils_solver = new ILS_CG(*gparams, ils_tol, ils_min_iter, ils_sigma, ils_max_iter, ils_tik_lam);
        }
        else if (gparams->ils_method == NLCG)
        {
            ils_solver = new ILS_NLCG(*gparams, ils_sigma, ils_max_iter, ils_tik_lam);
        }
        else if (gparams->ils_method == BiCGstabl)
        {
            ils_solver = new ILS_BiCGstabl(*gparams, ils_tol, ils_sigma, ils_max_iter, ils_tik_lam);
        }
        else
        {
            spdlog::error("SolverOSQP::solve()  Unknown Indirect Linear Solver method");
            return SolveResult{};
        }
        ils_solver->set_workspace(ws);

        total_Ax_size = 0;
        for (int i = 0; i < gparams->all_op.size(); i++)
        {
            total_Ax_size += gparams->all_op[i]->Ax_size;
        }
        r_primal.setZero(total_Ax_size);

        int total_feval = 0;
        for (iiter = 0; iiter < max_iter; ++iiter)
        {
            spdlog::trace("Starting OSQP iteration {:d} SolverOSQP::solve", iiter);

            if (iiter > 0)
            {
                Xhat = ils_solver->solve(X);
            }
            else
            {
                Xhat = X;
            }

            if (Xhat.array().isNaN().any())
            {
                spdlog::error("NaN detected in Xhat at iteration {:d}. Stopping solver.", iiter);
                break;
            };

            // Update all constraints (do prox operations)
            update(Xhat);

            X = gamma_x * Xhat + (1 - gamma_x) * X;

            get_residuals(X);

            if (extra_debug)
            {
                hist_X.push_back(X);
            }

            if ((logger(X) > 0) && (iiter > min_iter))
            {
                break;
            }

            total_feval += ils_solver->hist_n_iter.back();
            if (total_feval > max_feval)
            {
                spdlog::info("Maximum function evaluations reached");
                break;
            }
        }

        SolveResult result;
        result.X = X;
        result.n_iter = iiter;
        final_log(X, result);

        delete ils_solver;
        spdlog::trace("Finished SolverOSQP::solve");

        last_result = result;
        return result;
    }

    void SolverOSQP::update(Eigen::VectorXd &X)
    {
        spdlog::trace("Starting SolverOSQP::update");

        for (int i = 0; i < gparams->all_op.size(); i++)
        {
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
        }

        spdlog::trace("Finished SolverOSQP::update");
    }

    void SolverOSQP::get_residuals(Eigen::VectorXd &X)
    {
        // Update feasibility metrics
        for (int i = 0; i < gparams->all_op.size(); i++)
        {
            Operator *op = gparams->all_op[i].get();
            op->forward_op(X, op->Ax_temp);
            op->get_feas(op->Ax_temp);
            op->check(op->Ax_temp);
        }
    }

    void SolverOSQP::get_result_X(double **out, int &out_size)
    {
        out_size = last_result.X.size();
        *out = last_result.X.data();
    }

} // namespace Gropt
