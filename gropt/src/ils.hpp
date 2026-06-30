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
        std::vector<int> hist_n_iter;
        std::vector<double> hist_rnorm0;  // Initial residual norm ||b - A x0||  (warm-start residual)
        std::vector<double> hist_rnorm;  // final residual norm at the stopping point 
        std::vector<double> hist_bnorm0;  // ||b||  (RHS norm)
        std::vector<int> hist_neg_curv;  // 1 if a non-positive curvature direction (pAp <= 0) was hit
        // min/max Rayleigh quotient (pAp/pp) over the solve — this turned out to be of little use
        std::vector<double> hist_min_curv;
        std::vector<double> hist_max_curv;

        bool collect_debug = false;  // Set by the Solver = extra_debug.

        GroptParams *gparams;
        std::vector<const WorkspaceSolver*> ws;

        int n_iter;
        double sigma;
        double tik_lam;
        // Absolute residual floor (relative to ||b||) for the inner stop: keeps the warm-start-
        // relative criterion (tol*rnorm0) from chasing the residual down too much when close.
        double tol_abs_rel = 1e-8;

        IndirectLinearSolver(GroptParams &gparams, int _n_iter, double _sigma, double _tik_lam);
        ~IndirectLinearSolver() = default;

        void set_workspace(const std::vector<WorkspaceSolver*>& _ws) {
            ws.assign(_ws.begin(), _ws.end());
        }

        virtual Eigen::VectorXd solve(Eigen::VectorXd &x0);
        virtual void get_lhs(Eigen::VectorXd &x, Eigen::VectorXd &out);
        virtual void get_rhs(Eigen::VectorXd &x0, Eigen::VectorXd &out);
};

} // namespace Gropt

#endif
