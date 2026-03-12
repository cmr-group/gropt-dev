#ifndef EQUILIBRATE_H
#define EQUILIBRATE_H

#include "Eigen/Dense"

#include "gropt_params.hpp"

namespace Gropt {

class Operator; // Forward declaration

enum class NormType { L2, Inf };

void estimate_row_col_norms(GroptParams &gparams, int n_reps, NormType norm_type, Eigen::VectorXd &row_norms,
                            Eigen::VectorXd &col_norms);

void get_eq_vecs(GroptParams &gparams, Eigen::VectorXd &row_norms, Eigen::VectorXd &col_norms);

void rescale_eq_vecs(GroptParams &gparams, double row_scale, double col_scale);

void equilibrate(GroptParams &gparams, int n_iter, int n_reps);

double estimate_spec_norm(GroptParams &gparams, int n_iters);

double estimate_individual_spec_norm(GroptParams &gparams, int n_iters, int op_idx);

} // namespace Gropt

#endif
