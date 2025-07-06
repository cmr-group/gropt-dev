#ifndef SOLVER_H
#define SOLVER_H

#include <iostream> 
#include <string>
#include <numeric>
#include <vector>
#include "Eigen/Dense"

#include "gropt_params.hpp"
#include "ils.hpp"

namespace Gropt {

class Solver
{
    public:

        GroptParams *gparams;
        IndirectLinearSolver *ils_solver;

        int N_iter = 2000;
        int log_interval = 20;
        int min_iter = 0;

        double ils_tol = 1e-3;
        int ils_max_iter = 10;
        int ils_min_iter = 2;
        double ils_sigma = 1e-4;
        double ils_tik_lam = 1e-4;

        bool extra_debug = false;
        std::vector<Eigen::VectorXd> hist_X;
        std::vector<int> hist_cg_iter; 

        Solver(GroptParams &gparams);
        ~Solver() = default;

        virtual void solve();
        virtual int logger(Eigen::VectorXd &X);
        virtual void final_log(Eigen::VectorXd &X);
};

} // namespace Gropt

#endif