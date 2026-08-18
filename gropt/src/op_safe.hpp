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

    // Softabs smoothing of the |.| nonlinearity (slew units). 0 = hard abs / sign=+-1 (original). When
    // >0, |v| -> sqrt(v^2+eps^2) and sign -> v/sqrt(v^2+eps^2): the sign becomes a smooth ramp through 0
    // instead of a +-1 jump, so the frozen linearization the solver uses changes continuously instead of
    // flipping when a (filtered) slew sample crosses zero.
    double safe_eps = 0.0;

    // When true (set only for the duration of the inner CG via freeze_linearization), forward() does NOT
    // recapture the |.| signs -- it uses the frozen signs1/2/3 and applies them LINEARLY (sign*value
    // instead of abs), so the operator the CG sees is a fixed linear map. Restore to false afterward.
    bool freeze_signs = false;

    Op_SAFE(const ProblemData &_pdata, double _stim_thresh, double _weight_mod);
    Op_SAFE(const ProblemData &_pdata, const Eigen::VectorXd &_stim_thresh_vec, double _weight_mod);
    virtual void init();

    virtual void forward(Eigen::VectorXd &X, Eigen::VectorXd &out);
    virtual void transpose(Eigen::VectorXd &X, Eigen::VectorXd &out);
    virtual void prox(Eigen::VectorXd &X);
    virtual void check(Eigen::VectorXd &X);

    // Freeze the |.| linearization at X for the inner CG (recapture signs once, then hold), and restore.
    virtual void freeze_linearization(Eigen::VectorXd &X) override;
    virtual void unfreeze_linearization() override { freeze_signs = false; }
    virtual double linearization_error(const Eigen::VectorXd &x_new) override;
    virtual double constraint_violation(const Eigen::VectorXd &x_new) override;

    // SAFE stacks n_terms blocks per axis (Ax_size = n_terms*Naxis*N); each N-block is a separate
    // time series, so it must be resized block-by-block rather than as one Ax_size/Naxis chunk.
    virtual std::vector<int> Ax_block_lengths() const override { return std::vector<int>(n_terms * Naxis, N); }
};

} // namespace Gropt

#endif
