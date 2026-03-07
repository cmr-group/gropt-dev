#ifndef GROPT_UTILS_H
#define GROPT_UTILS_H

/**
 * In general this is just a place for random usages of the GrOpt
 * operators for applications that aren't actual optimization.
 *
 * i.e. GIRF respone calculation or spectral calcs, or getting a
 * PNS curve for a waveform
 */

#include "Eigen/Dense"
#include <iostream>
#include <string>
#include <vector>

namespace Gropt {

// Eigen overload: demo params
Eigen::VectorXd get_SAFE_eigen(const Eigen::VectorXd &G, int Naxis, double dt, bool true_safe, int new_first_axis);
// Eigen overload: custom params
Eigen::VectorXd get_SAFE_eigen(const Eigen::VectorXd &G, int Naxis, double dt, bool true_safe, int new_first_axis,
                               const Eigen::VectorXd &tau1, const Eigen::VectorXd &tau2, const Eigen::VectorXd &tau3,
                               const Eigen::VectorXd &a1, const Eigen::VectorXd &a2, const Eigen::VectorXd &a3,
                               const Eigen::VectorXd &stim_limit, const Eigen::VectorXd &g_scale);

void test_eigen_assertions(int test_type);

} // namespace Gropt

#endif