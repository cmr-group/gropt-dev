#ifndef SOLVER_H
#define SOLVER_H

#include "Eigen/Dense"
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include "gropt_params.hpp"
#include "ils.hpp"
#include "workspace_solver.hpp"

namespace Gropt {

class GroptParams; // Forward declaration of GroptParams class

struct DebugSolver {
    std::vector<Eigen::VectorXd> hist_X;
    std::vector<Eigen::VectorXd> hist_Ax;
    std::vector<Eigen::VectorXd> hist_z;
    std::vector<Eigen::VectorXd> hist_y;
    std::vector<Eigen::VectorXd> hist_Aty;
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
    int min_iter = 0;
    double gamma_x = 1.6;

    double ils_tol = 1e-3;
    int ils_max_iter = 10;
    int ils_min_iter = 2;
    double ils_sigma = 1e-4;
    double ils_tik_lam = 1e-4;

    bool extra_debug = false;
    DebugSolver debug_solver;

    std::vector<int> hist_cg_iter;

    // Iteration counter (was on GroptParams)
    int iiter = 0;

    Solver() = default;
    ~Solver() = default;

    virtual SolveResult solve(GroptParams &_gparams);
    virtual int logger(Eigen::VectorXd &X);
    virtual void final_log(Eigen::VectorXd &X, SolveResult &result);
    virtual void set_general_params(int min_iter, int max_iter, int log_interval, double gamma_x, int max_feval);
    virtual void set_ils_params(double ils_tol, int ils_max_iter, int ils_min_iter, double ils_sigma,
                                double ils_tik_lam);
};

} // namespace Gropt

#endif
