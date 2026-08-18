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

        int rw_interval = 8;
        double rw_e_corr = 0.4;
        double rw_eps = 1e-36;
        double rw_scalelim = 1.5;

        int grw_min_infeasible = 20;
        int grw_interval = 20;
        double grw_mod = 2.0;
        
        // true = grw keeps constraint geometric mean constant, false = constantly raises the worst constraint weight (default). 
        // The latter is more aggressive, but can drown the objective if the constraints are many and strong.
        bool grw_balanced = false;

        // Re-project the over-relaxed iterate onto the equality (moment/eddy/concomitant) surface every
        // outer iteration.
        bool reproject_iterate = true;

        // Low-frequency projection (fft_tools). Each outer iteration, project the iterate onto
        // frequencies <= cutoff_freq [Hz] (per-axis DCT hard low-pass) to suppress high-frequency
        // oscillation.
        double cutoff_freq = -1.0;
        int cutoff_iter = -1;
        // Raised-cosine roll-off width of the low-pass, as a fraction of the cutoff bin. 0 = hard cutoff
        double cutoff_trans = 0.0;

        // Trust-region step control (default off). After the inner CG, a StepMonitor tests whether the
        // linearized model held over the step; if not, re-solve from X with the proximal sigma scaled up.
        bool tr_enable = false;
        double tr_tol = -1.0;     // monitor reject threshold; <=0 keeps the monitor's own default
                                  // (0.2 linearization_error, 0.02 feasibility) since the good tol differs
        double tr_bump = 4.0;     // proximal-sigma multiplier per reject
        int tr_max_reject = 5;    // max re-solves per outer iteration before taking the damped step anyway
        double tr_decay = 0.5;    // sigma relaxation toward ils_sigma on an accepted step
        std::string tr_monitor = "linearization_error";
        std::unique_ptr<StepMonitor> step_monitor; // built in solve() when tr_enable

        // Feasibility-gated objective (default off). Each outer iteration the objective pull is scaled by
        // exp(-total_constraint_violation / obj_gate_scale): ~0 while infeasible (no cold overshoot past
        // the constraints), ->1 when feasible (climb to the max b). obj_gate_scale sets how sharply the
        // gate opens as the violation shrinks (in constraint units, e.g. ~0.05 of the SAFE limit).
        bool obj_gate_enable = false;
        double obj_gate_scale = 0.05;

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
