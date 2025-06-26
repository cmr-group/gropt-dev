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

        SolverGroptSDMM(GroptParams &gparams) : Solver(gparams){gparams.solver_method = GROPT_SDMM;};

        double gamma_x = 1.6;
        
        int total_Ax_size;

        Eigen::VectorXd Px;
        Eigen::VectorXd r_dual; 
        Eigen::VectorXd r_primal;  

        virtual void solve();
        void update(Eigen::VectorXd &X);
        void get_residuals(Eigen::VectorXd &X);
};

} // namespace Gropt

#endif