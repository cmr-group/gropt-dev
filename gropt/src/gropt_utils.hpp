#ifndef GROPT_UTILS_H
#define GROPT_UTILS_H

/**
 * Non-optimization uses of the GrOpt operators, e.g. the SAFE PNS curve of a waveform.
 */

#include "Eigen/Dense"
#include <iostream>
#include <string>
#include <vector>

namespace Gropt {

// SAFE PNS curve of G (axis-major), demo params
Eigen::VectorXd get_SAFE_eigen(const Eigen::VectorXd &G, int Naxis, double dt, bool true_safe, int new_first_axis);
// Same, custom params
Eigen::VectorXd get_SAFE_eigen(const Eigen::VectorXd &G, int Naxis, double dt, bool true_safe, int new_first_axis,
                               const Eigen::VectorXd &tau1, const Eigen::VectorXd &tau2, const Eigen::VectorXd &tau3,
                               const Eigen::VectorXd &a1, const Eigen::VectorXd &a2, const Eigen::VectorXd &a3,
                               const Eigen::VectorXd &stim_limit, const Eigen::VectorXd &g_scale);

void test_eigen_assertions(int test_type);

} // namespace Gropt

#endif