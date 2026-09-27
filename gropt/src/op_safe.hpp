#ifndef OP_SAFE_H
#define OP_SAFE_H

/**
 * Constraint on the SAFE-model PNS prediction of the gradient waveform.
 *
 * Per sample the axes combine as the root-sum-square, sqrt(sum_j (stim_j / thresh_j)^2) <= 1, as the scanner
 * adds them (per-axis limits would be up to sqrt(Naxis) weaker). Not rotationally invariant (per-axis
 * constants differ), so rot_variant must be true.
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

    // One-pole filter step. false: alpha = dt/(tau+dt), backward Euler as in the reference Matlab; raster
    // dependent (~8% low at 400 us). true: 1 - exp(-dt/tau), exact for piecewise-linear waveforms and raster
    // invariant (conservative; ~0.2% apart at 10 us).
    bool alpha_exact = false;

    // Filter coefficients, derived by calc_alphas(dt) rather than set by the user
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
    bool thresh_from_vec = false; // a caller-supplied limit vector: must be Naxis*N long, never rebuilt

    Eigen::VectorXd signs1;
    Eigen::VectorXd signs2;
    Eigen::VectorXd signs3;
    Eigen::VectorXd slew_temp;

    // One-pole IIR pieces, shared with Op_SAFE_Slack so the two forwards differ only in term 2.
    void compute_slew(const Eigen::VectorXd &X);                        // -> slew_temp
    void lowpass(Eigen::VectorXd &v, const std::vector<double> &alpha) const;   // causal, in place
    void lowpass_T(Eigen::VectorXd &v, const std::vector<double> &alpha) const; // its adjoint
    void take_abs(Eigen::VectorXd &v, Eigen::VectorXd &signs) const;    // |.|, softabs or frozen sign
    void diff_T(Eigen::VectorXd &out) const;                            // D^T on the waveform block
    void emit(Eigen::VectorXd &out) const;   // out = [a1 stim1; a2 stim2; a3 stim3] * g_scale / stim_limit
    void absorb(const Eigen::VectorXd &X);   // its adjoint, into stim1/2/3

    Eigen::VectorXd stim1;
    Eigen::VectorXd stim2;
    Eigen::VectorXd stim3;

  public:
    SAFEParams safe_params;
    int n_terms = 3;

    // Softabs smoothing of |.| [T/m/s]: |v| -> sqrt(v^2 + eps^2); 0 = exact |.|
    double safe_eps = 0.0;

    // Emit terms 1/3 signed; prox() takes their |.| exactly (abs outside the filter, no sign freezing).
    // Term 2's abs is inside the filter, so it needs Op_SAFE_Slack instead.
    bool signed_terms13 = false;

    // Set during the inner CG: forward() applies the held signs1/2/3 linearly instead of |.|
    bool freeze_signs = false;

    // Summed stimulation of axis j at sample i, from the n_terms blocks of X (Ax space)
    double axis_stim(const Eigen::VectorXd &X, int j, int i) const;
    // Same, for scoring (check, constraint_violation): |.| on signed terms 1/3 so they cannot cancel term 2
    double stim_mag(const Eigen::VectorXd &X, int j, int i) const;

    Op_SAFE(const ProblemData &_pdata, double _stim_thresh, double _weight_mod);
    Op_SAFE(const ProblemData &_pdata, const Eigen::VectorXd &_stim_thresh_vec, double _weight_mod);
    virtual void init();

    virtual void forward(Eigen::VectorXd &X, Eigen::VectorXd &out);
    virtual void transpose(Eigen::VectorXd &X, Eigen::VectorXd &out);
    virtual void prox(Eigen::VectorXd &X);
    void prox_signed(Eigen::VectorXd &X); // projection in magnitude space, used when signed_terms13

    virtual void check(Eigen::VectorXd &X);

    // Hold the |.| signs captured at X during the inner CG
    virtual void freeze_linearization(Eigen::VectorXd &X) override;
    virtual void unfreeze_linearization() override { freeze_signs = false; }
    virtual double linearization_error(const Eigen::VectorXd &x_new) override;
    virtual double constraint_violation(const Eigen::VectorXd &x_new) override;

    // Ax is n_terms*Naxis separate length-N time series, not Naxis blocks
    virtual std::vector<int> Ax_block_lengths() const override { return std::vector<int>(n_terms * Naxis, N); }
};

} // namespace Gropt

#endif
