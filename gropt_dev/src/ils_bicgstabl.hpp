#ifndef ILS_BICGSTABL_H
#define ILS_BICGSTABL_H

#include <iostream> 
#include <string>
#include <vector>
#include <cmath>
#include "Eigen/Dense"

#include "ils.hpp"

namespace Gropt {

class ILS_BiCGstabl : public IndirectLinearSolver
{
    public:
        int ell = 2;

        Eigen::MatrixXd rs;
        Eigen::MatrixXd us;
        Eigen::MatrixXd tau;

        Eigen::VectorXd x;
        Eigen::VectorXd r_shadow;
        Eigen::VectorXd r;
        
        Eigen::VectorXd Ax;
        Eigen::VectorXd b;
        
        Eigen::VectorXd sigma;
        Eigen::VectorXd gamma;  
        Eigen::VectorXd gammap;
        Eigen::VectorXd gammapp;

        double tol;

        ILS_BiCGstabl(GroptParams &_gparams, double _tol, double _sigma, int _n_iter, double _tik_lam);

        // Runs conventional conjugate gradient
        Eigen::VectorXd solve(Eigen::VectorXd &x_in) override;

};

}  // end namespace Gropt

#endif