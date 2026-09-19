#ifndef ILS_H
#define ILS_H

#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include "Eigen/Dense"

#include "gropt_params.hpp"
#include "workspace_solver.hpp"

namespace Gropt {

class IndirectLinearSolver
{
    public:
        std::string name;
        std::chrono::steady_clock::time_point start_time;
        std::chrono::steady_clock::time_point stop_time;
        std::chrono::duration<double, std::micro> elapsed_us;
        // Per-solve histories; only hist_n_iter has a leading -1 placeholder.
        std::vector<int> hist_n_iter;
        std::vector<double> hist_rnorm0;  // initial residual ||b - A x0||
        std::vector<double> hist_rnorm;  // final residual
        std::vector<double> hist_bnorm0;  // ||b||

        GroptParams *gparams;
        std::vector<const WorkspaceSolver*> ws;

        int n_iter;
        double sigma;
        double tik_lam;
        // Stop-threshold floor relative to ||b||: CG/BiCGstabl stop at max(tol*rnorm0, tol_abs_rel*||b||),
        // so a tiny warm-start residual cannot set an unreachable target.
        double tol_abs_rel = 1e-8;

        IndirectLinearSolver(GroptParams &gparams, int _n_iter, double _sigma, double _tik_lam);
        // Virtual: ILS objects are deleted through the base pointer Solver::ils_solver.
        virtual ~IndirectLinearSolver() = default;

        void set_workspace(const std::vector<WorkspaceSolver*>& _ws) {
            ws.assign(_ws.begin(), _ws.end());
        }

        virtual Eigen::VectorXd solve(Eigen::VectorXd &x0);
        virtual void get_lhs(Eigen::VectorXd &x, Eigen::VectorXd &out);
        virtual void get_rhs(Eigen::VectorXd &x0, Eigen::VectorXd &out);
};

} // namespace Gropt

#endif
