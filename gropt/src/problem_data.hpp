#ifndef PROBLEM_DATA_H
#define PROBLEM_DATA_H

#include "Eigen/Dense"

namespace Gropt {

struct ProblemData {
    int N = 10;
    int Naxis = 1;
    double dt = 10e-6;

    Eigen::VectorXd X0;       // Initial guess
    Eigen::VectorXd inv_vec;  // Inversion vector (diffusion encoding sign flips)
    Eigen::VectorXd set_vals; // Fixed values (NaN = free)
    Eigen::VectorXd fixer;    // Binary mask: 0 = fixed, 1 = free
};

}  // namespace Gropt

#endif
