#ifndef ILS_NLCG_H
#define ILS_NLCG_H

#include <iostream> 
#include <string>
#include <vector>
#include <cmath>
#include "Eigen/Dense"

#include "ils.hpp"

namespace Gropt {

class ILS_NLCG : public IndirectLinearSolver
{
    public:
        Eigen::VectorXd x0;
        Eigen::VectorXd x1;
        Eigen::VectorXd b;
        Eigen::VectorXd Ax;
        Eigen::VectorXd r;

        double eta = 0.8;
        double theta = 0.5;


        ILS_NLCG(GroptParams &_gparams, double _sigma, int _n_iter, double _tik_lam);

        // Runs conventional conjugate gradient
        Eigen::VectorXd solve(Eigen::VectorXd &x0) override;

};

}  // end namespace Gropt

#endif