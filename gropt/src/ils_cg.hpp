#ifndef ILS_CG_H
#define ILS_CG_H

#include <iostream> 
#include <string>
#include <vector>
#include "Eigen/Dense"

#include "ils.hpp"

namespace Gropt {

class ILS_CG : public IndirectLinearSolver
{
    public:
        Eigen::VectorXd x;
        Eigen::VectorXd b;
        Eigen::VectorXd Ax;
        Eigen::VectorXd Ap;
        Eigen::VectorXd r;
        Eigen::VectorXd p;
        
        double tol;
        int min_iter;

        ILS_CG(GroptParams &_gparams, double _tol, int _min_iter, double _sigma, int _n_iter, double _tik_lam);

        // Runs conventional conjugate gradient
        Eigen::VectorXd solve(Eigen::VectorXd &x0) override;

};

}  // end namespace Gropt

#endif