#ifndef ILS_H
#define ILS_H

#include <iostream> 
#include <string>
#include <vector>
#include <chrono>
#include "Eigen/Dense"

#include "gropt_params.hpp"

namespace Gropt {

class IndirectLinearSolver
{
    public:
        std::string name;
        std::chrono::steady_clock::time_point start_time;
        std::chrono::steady_clock::time_point stop_time;
        std::chrono::duration<double, std::micro> elapsed_us;
        std::vector<int> hist_n_iter; 

        GroptParams *gparams;

        int n_iter;
        double sigma;
        double tik_lam;

        IndirectLinearSolver(GroptParams &gparams, int _n_iter, double _sigma, double _tik_lam);
        ~IndirectLinearSolver() = default;

        virtual Eigen::VectorXd solve(Eigen::VectorXd &x0);
        virtual Eigen::VectorXd get_lhs(Eigen::VectorXd &x, Eigen::VectorXd &out);
        virtual Eigen::VectorXd get_rhs(Eigen::VectorXd &x0, Eigen::VectorXd &out);
};

} // namespace Gropt

#endif