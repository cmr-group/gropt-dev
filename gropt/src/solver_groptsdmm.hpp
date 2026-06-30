#ifndef SOLVER_GROPTSDMM_H
#define SOLVER_GROPTSDMM_H

#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include "Eigen/Dense"

#include "solver.hpp"
#include "gropt_params.hpp"
#include "workspace_sdmm.hpp"

namespace Gropt
{

    class Op_BValue; // forward decl (record_debug reads the achieved b-value)

    class SolverGroptSDMM : public Solver
    {
    public:
        SolverGroptSDMM() = default;

        // Per-operator SDMM workspaces (typed; base class Solver::ws points into this)
        std::vector<WorkspaceSDMM> sdmm_ws;

        int total_Ax_size;

        int rw_interval = 8;
        double rw_e_corr = 0.4;
        double rw_eps = 1e-36;
        double rw_scalelim = 1.5;

        int grw_min_infeasible = 20;
        int grw_interval = 20;
        double grw_mod = 2.0;

        Eigen::VectorXd Px;
        Eigen::VectorXd r_dual;
        Eigen::VectorXd r_primal;

        virtual SolveResult solve(GroptParams &_gparams);
        void update(Eigen::VectorXd &X);
        void get_residuals(Eigen::VectorXd &X);
        void set_sdmm_params(int rw_interval, double rw_e_corr, double rw_eps, double rw_scalelim,
                             int grw_min_infeasible, int grw_interval, double grw_mod);

      private:
        // solve() helpers -- behavior-preserving extractions that keep solve() readable.
        Eigen::VectorXd resolve_initial_primal();             // warm-start resize, or the cold X0
        void init_workspaces(const Eigen::VectorXd &X0_init); // per-op SDMM workspaces + warm-dual injection
        void record_debug(Eigen::VectorXd &X, Op_BValue *bval_op); // per-iteration extra_debug capture
    };

} // namespace Gropt

#endif
