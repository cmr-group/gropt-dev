#ifndef OP_SAFE_H
#define OP_SAFE_H

/**
 * Constriant on SAFE model prediction of waveforms.
 */

#include "Eigen/Dense"
#include <iostream>
#include <math.h>
#include <string>

#include "op_main.hpp"

namespace Gropt {

class SAFEParams {

  public:
    std::vector<double> tau1 = std::vector<double>(3, 0.0);
    std::vector<double> tau2 = std::vector<double>(3, 0.0);
    std::vector<double> tau3 = std::vector<double>(3, 0.0);
    std::vector<double> a1 = std::vector<double>(3, 0.0);
    std::vector<double> a2 = std::vector<double>(3, 0.0);
    std::vector<double> a3 = std::vector<double>(3, 0.0);
    std::vector<double> stim_limit = std::vector<double>(3, 0.0);
    std::vector<double> g_scale = std::vector<double>(3, 0.0);

    // These aren't real parameters because they depend on dt, maybe they should be in Op_SAFE?
    std::vector<double> alpha1 = std::vector<double>(3, 0.0);
    std::vector<double> alpha2 = std::vector<double>(3, 0.0);
    std::vector<double> alpha3 = std::vector<double>(3, 0.0);

    SAFEParams() = default;
    void set_demo_params();
    void set_params(const Eigen::VectorXd &_tau1, const Eigen::VectorXd &_tau2, const Eigen::VectorXd &_tau3,
                    const Eigen::VectorXd &_a1, const Eigen::VectorXd &_a2, const Eigen::VectorXd &_a3,
                    const Eigen::VectorXd &_stim_limit, const Eigen::VectorXd &_g_scale);
    void calc_alphas(double dt);
    void swap_first_axes(int new_first_axis);
};

class Op_SAFE : public Operator {
  protected:
    double stim_thresh;

    Eigen::VectorXd stim_thresh_vec;

    Eigen::VectorXd signs1;
    Eigen::VectorXd signs2;
    Eigen::VectorXd signs3;
    Eigen::VectorXd stim1;
    Eigen::VectorXd stim2;
    Eigen::VectorXd stim3;

  public:
    SAFEParams safe_params;
    int n_terms = 3;

    Op_SAFE(const ProblemData &_pdata, double _stim_thresh, double _weight_mod);
    Op_SAFE(const ProblemData &_pdata, const Eigen::VectorXd &_stim_thresh_vec, double _weight_mod);
    virtual void init();

    virtual void forward(Eigen::VectorXd &X, Eigen::VectorXd &out);
    virtual void transpose(Eigen::VectorXd &X, Eigen::VectorXd &out);
    virtual void prox(Eigen::VectorXd &X);
    virtual void check(Eigen::VectorXd &X);

    // SAFE stacks n_terms blocks per axis (Ax_size = n_terms*Naxis*N); each N-block is a separate
    // time series, so it must be resized block-by-block rather than as one Ax_size/Naxis chunk.
    virtual std::vector<int> Ax_block_lengths() const override { return std::vector<int>(n_terms * Naxis, N); }
};

} // namespace Gropt

#endif
