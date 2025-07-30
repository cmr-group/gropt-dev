#ifndef SOLVER_GROPTSDMM_H
#define SOLVER_GROPTSDMM_H

#include <iostream> 
#include <algorithm> 
#include <string>
#include <vector>
#include "Eigen/Dense"

#include "solver.hpp"
#include "gropt_params.hpp"

namespace Gropt {

class SolverGroptSDMM : public Solver
{
    public:

        SolverGroptSDMM() = default;

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

        virtual void solve(GroptParams &_gparams);
        void update(Eigen::VectorXd &X);
        void get_residuals(Eigen::VectorXd &X);
        void set_sdmm_params(int rw_interval, double rw_e_corr, double rw_eps, double rw_scalelim,
                             int grw_min_infeasible, int grw_interval, double grw_mod);
};

} // namespace Gropt

#endif