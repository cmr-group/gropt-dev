#ifndef GROPT_UTILS_H
#define GROPT_UTILS_H

/**
 * In general this is just a place for random usages of the GrOpt
 * operators for applications that aren't actual optimization.
 * 
 * i.e. GIRF respone calculation or spectral calcs, or getting a 
 * PNS curve for a waveform
 */

#include <iostream> 
#include <string>
#include <vector>
#include "Eigen/Dense"

namespace Gropt {
    
void get_SAFE(int N, int Naxis, double dt, double *G_in,
              bool true_safe, int new_first_axis, bool demo_params,
              double *tau1, double *tau2, double *tau3,
              double *a1, double *a2, double *a3,
              double *stim_limit, double *g_scale,
              double **out, int &out_size);

// Eigen overload: demo params
Eigen::VectorXd get_SAFE_eigen(const Eigen::VectorXd &G, int Naxis, double dt,
                               bool true_safe, int new_first_axis);
// Eigen overload: custom params
Eigen::VectorXd get_SAFE_eigen(const Eigen::VectorXd &G, int Naxis, double dt,
                               bool true_safe, int new_first_axis,
                               const Eigen::VectorXd &tau1, const Eigen::VectorXd &tau2, const Eigen::VectorXd &tau3,
                               const Eigen::VectorXd &a1, const Eigen::VectorXd &a2, const Eigen::VectorXd &a3,
                               const Eigen::VectorXd &stim_limit, const Eigen::VectorXd &g_scale);

}  // close "namespace Gropt"

#endif