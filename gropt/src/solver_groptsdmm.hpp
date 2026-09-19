#ifndef SOLVER_GROPTSDMM_H
#define SOLVER_GROPTSDMM_H

#include <iostream>
#include <algorithm>
#include <memory>
#include <string>
#include <vector>
#include "Eigen/Dense"

#include "solver.hpp"
#include "gropt_params.hpp"
#include "step_monitor.hpp"
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

        bool bb_enable = true;   // per-operator BB (spectral) reweighting every rw_interval iterations
        bool grw_enable = true;  // global reweighting of the worst persistently infeasible operator

        int rw_interval = 8;
        double rw_e_corr = 0.4;
        double rw_eps = 1e-36;
        double rw_scalelim = 1.5;

        int grw_min_infeasible = 20;
        int grw_interval = 20;
        double grw_mod = 2.0;
        // true: rescale all constraint weights after each bump to keep their geometric mean fixed.
        // false (default): only raise the worst one; more aggressive, can drown out the objective.
        bool grw_balanced = false;

        // Re-project the over-relaxed iterate onto the equality surface every outer iteration.
        bool reproject_iterate = true;

        // Low-frequency projection (fft_tools): each outer iteration, low-pass the iterate at cutoff_freq [Hz]
        // with a per-free-run DST-I to suppress high-frequency oscillation. <= 0 disables.
        double cutoff_freq = -1.0;
        int cutoff_iter = -1;
        // Raised-cosine roll-off width as a fraction of the cutoff bin; 0 = brick wall. Wider suppresses more
        // near-cutoff oscillation but strips harmonics that keep plateaus flat, so ripple can go either way.
        double cutoff_trans = 0.0;

        // Trust-region step control (default off): when the StepMonitor rejects an inner-CG step, re-solve
        // from X with the proximal sigma scaled up.
        bool tr_enable = false;
        double tr_tol = -1.0;     // monitor reject threshold; <= 0 uses the monitor's default
        double tr_bump = 4.0;     // proximal-sigma multiplier per reject
        int tr_max_reject = 5;    // max re-solves per outer iteration before taking the damped step anyway
        double tr_decay = 0.5;    // sigma relaxation toward ils_sigma on an accepted step
        std::string tr_monitor = "linearization_error";

        // Feasibility-gated objective (default off): each outer iteration the objective pull is scaled by
        // exp(-violation / obj_gate_scale), violation = sum of Operator::constraint_violation over all_op.
        // Only Op_SAFE implements constraint_violation, so the gate currently reacts to PNS/CNS only.
        bool obj_gate_enable = false;
        double obj_gate_scale = 0.05;

        virtual SolveResult solve(GroptParams &_gparams);
        void update(Eigen::VectorXd &X);
        void get_residuals(Eigen::VectorXd &X);
        void set_sdmm_params(int rw_interval, double rw_e_corr, double rw_eps, double rw_scalelim,
                             int grw_min_infeasible, int grw_interval, double grw_mod);

      private:
        // solve() helpers
        Eigen::VectorXd resolve_initial_primal();             // warm-start resize, or the cold X0
        void init_workspaces(const Eigen::VectorXd &X0_init); // per-op SDMM workspaces + warm-dual injection
        void record_debug(Eigen::VectorXd &X, Op_BValue *bval_op); // per-iteration extra_debug capture
    };

} // namespace Gropt

#endif
