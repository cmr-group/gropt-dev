#ifndef OP_SAFE_SLACK_H
#define OP_SAFE_SLACK_H

/**
 * SAFE with its abs-inside-the-filter term lifted: term 2 reads an auxiliary u >= |slew| (kept honest by
 * Op_SlackAbs) instead of computing LP2(|slew|) itself, which makes that term linear.
 *
 * Only term 2, because only term 2 needs it. Op_SAFE freezes the sign of each |.| for the inner CG, so
 * the operator the CG solves with differs from the one the z/y update applies; measured per term on the
 * diffusion blow-ups that mismatch is 1.30 for term 2 against 0.000 for term 1, since the low-pass makes
 * the FILTERED slew's sign stable across a step while the RAW slew's flips on 11-36% of samples. Terms 1
 * and 3 have their abs outside the filter, so Op_SAFE::signed_terms13 handles them without a variable.
 *
 * A separate class so Op_SAFE, the path every existing problem uses, is untouched. Ax keeps the same
 * layout, so prox() and the combined root-sum-square limit are inherited unchanged.
 */

#include <string>

#include "op_safe.hpp"

namespace Gropt {

// Name of the shared u >= |slew| block. Every Op_SAFE_Slack reads the same one: the models differ
// in their filters and coefficients, but the slew they bound is the waveform's, not the model's.
inline const char *SAFE_SLACK_BLOCK = "abs_slew";

class Op_SAFE_Slack : public Op_SAFE {
  protected:
    std::string aux_name;
    int u_offset = -1;
    Eigen::VectorXd x_cache; // waveform of the last forward(), for scoring the true model in check()

  public:
    Op_SAFE_Slack(const ProblemData &_pdata, double _stim_thresh, double _weight_mod,
                  const std::string &_aux_name);
    Op_SAFE_Slack(const ProblemData &_pdata, const Eigen::VectorXd &_stim_thresh_vec, double _weight_mod,
                  const std::string &_aux_name);

    const std::string &aux_block() const { return aux_name; }

    // signed_terms13 lives on Op_SAFE: terms 1 and 3 are not what the lifting is about.

    virtual void declare_aux(ProblemData &pd) override;
    virtual void init() override;

    virtual void forward(Eigen::VectorXd &X, Eigen::VectorXd &out) override;
    virtual void transpose(Eigen::VectorXd &X, Eigen::VectorXd &out) override;
    virtual void check(Eigen::VectorXd &X) override;
    virtual double constraint_violation(const Eigen::VectorXd &x_new) override;

  protected:
    // True (unlifted) SAFE response of the cached waveform, in Ax layout and physical units.
    void true_Ax(Eigen::VectorXd &out);
};

} // namespace Gropt

#endif
