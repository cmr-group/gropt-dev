#ifndef SOLVER_H
#define SOLVER_H

#include "Eigen/Dense"
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include "gropt_params.hpp"
#include "ils.hpp"
#include "warmstart.hpp"
#include "workspace_solver.hpp"

namespace Gropt {

class GroptParams; // Forward declaration of GroptParams class

// Per-iteration solver traces, filled only when Solver::extra_debug is true.
struct DebugSolver {
    std::vector<Eigen::VectorXd> hist_X;
    std::vector<Eigen::VectorXd> hist_Ax;
    std::vector<Eigen::VectorXd> hist_z;
    std::vector<Eigen::VectorXd> hist_y;
    std::vector<Eigen::VectorXd> hist_Aty;
    std::vector<std::vector<double>> hist_weight;
    std::vector<std::vector<double>> hist_gamma;
    std::vector<double> hist_gamma_x;
    // Boyd ADMM residuals per iteration, per non-projected operator:
    //   primal = ||A x - z||                 (consensus gap)
    //   dual   = ||rho * A^T (z - z_prev)||  (stationarity)
    std::vector<std::vector<double>> hist_r_prim;
    std::vector<std::vector<double>> hist_r_dual;
    // Per-operator feasibility per iteration:
    //   r_feas = relative distance of A x from the feasible set
    //   feas   = binary feasible flag (1/0)
    std::vector<std::vector<double>> hist_r_feas;
    std::vector<std::vector<int>> hist_feas;
    // 1 if all operators were feasible this iteration, else 0.
    std::vector<int> hist_all_feas;
    // Achieved b-value per iteration. Empty if no b-value operator is present.
    std::vector<double> hist_bvalue;
    // Outer iteration of the returned (best feasible) iterate, or -1 if none was feasible.
    int best_feasible_iter = -1;
    // Inner linear-solver diagnostics, one entry per inner solve (hist_cg_iter[0] is a -1 placeholder):
    //   n_iter = iterations the inner solver took
    //   rnorm0 = initial residual ||b - A x0|| (warm-start residual)
    //   rnorm  = final residual at the inner-solver stop
    //   bnorm0 = ||b|| (RHS norm)
    std::vector<int> hist_cg_iter;
    std::vector<double> hist_cg_rnorm0;
    std::vector<double> hist_cg_rnorm;
    std::vector<double> hist_cg_bnorm0;
    // Objective vs constraint pull on the x-update RHS, per outer iteration:
    std::vector<double> hist_obj_pull;                 // ||g_obj||  linearized objective pull
    std::vector<double> hist_con_pull;                 // ||Σ Aᵀy||  total constraint pull
    std::vector<std::vector<double>> hist_con_pull_op; // per-operator ||Aᵀy||
};

class Solver {
  public:
    GroptParams *gparams;
    IndirectLinearSolver *ils_solver;

    // Per-operator workspaces (pointers into subclass-owned typed vectors)
    std::vector<WorkspaceSolver *> ws;

    int max_iter = 2000;
    int max_feval = 12000;
    int log_interval = 20;
    int min_iter = 1;
    double gamma_x = 1.6;
    int obj_patience = 20;  // stop after this many feasible iters with no objective improvement
    double obj_rtol = 1e-4; // relative objective-improvement threshold

    double ils_tol = 1e-3;
    int ils_max_iter = 20;
    int ils_min_iter = 2;
    double ils_sigma = 1e-4;
    double ils_tik_lam = 0.0;

    bool extra_debug = false;
    DebugSolver debug_solver;

    // Warm-start state (see warmstart.hpp).
    WarmStart warmstart;       // input: loaded by set_warmstart(), consumed by solve()
    WarmStart best_warmstart;  // output: captured at the best feasible iterate, returned by get_warmstart()

    // Outer iteration counter
    int iiter = 0;

    Solver() = default;
    virtual ~Solver() = default; // deleted through Solver* when the solver type is chosen at runtime

    virtual SolveResult solve(GroptParams &_gparams);
    virtual int logger(Eigen::VectorXd &X);
    virtual void final_log(Eigen::VectorXd &X, SolveResult &result);

    // Build a warm-start snapshot from the current per-operator workspaces (ws) and primal X.
    WarmStart capture_warmstart(const Eigen::VectorXd &X);
    // Snapshot captured at the returned best-feasible iterate (or the final iterate if none feasible).
    WarmStart get_warmstart() { return best_warmstart; }
    // Load a snapshot to warm-start the next solve().
    void set_warmstart(const WarmStart &ws) {
        warmstart = ws;
        warmstart.active = true;
    }
    virtual void set_general_params(int min_iter, int max_iter, int log_interval, double gamma_x, int max_feval,
                                    int obj_patience);
    virtual void set_ils_params(double ils_tol, int ils_max_iter, int ils_min_iter, double ils_sigma,
                                double ils_tik_lam);
};

} // namespace Gropt

#endif
