#ifndef SOLVER_OSQP_H
#define SOLVER_OSQP_H

#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include "Eigen/Dense"

#include "solver.hpp"
#include "gropt_params.hpp"
#include "workspace_osqp.hpp"

namespace Gropt
{

    class SolverOSQP : public Solver
    {
    public:
        SolverOSQP() = default;

        // Per-operator OSQP workspaces (typed; base class Solver::ws points into this)
        std::vector<WorkspaceOSQP> osqp_ws;

        int total_Ax_size;

        Eigen::VectorXd Px;
        Eigen::VectorXd r_dual;
        Eigen::VectorXd r_primal;

        virtual SolveResult solve(GroptParams &_gparams);
        void update(Eigen::VectorXd &X);
        void get_residuals(Eigen::VectorXd &X);
    };

} // namespace Gropt

#endif
