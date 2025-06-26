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
        Eigen::VectorXd x;
        Eigen::VectorXd b;
        Eigen::VectorXd Ax;
        Eigen::VectorXd dx;
        Eigen::VectorXd s0;
        Eigen::VectorXd s1;
        Eigen::VectorXd r;

        double alpha;
        double beta;


        ILS_NLCG(GroptParams &_gparams);

        // Runs conventional conjugate gradient
        Eigen::VectorXd solve(Eigen::VectorXd &x0) override;
        double line_search();

};

}  // end namespace Gropt

#endif