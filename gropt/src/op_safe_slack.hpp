#ifndef OP_SAFE_SLACK_H
#define OP_SAFE_SLACK_H

/**
 * SAFE with its abs-inside-the-filter term lifted: term 2 reads an auxiliary u >= |slew| (kept honest by
 * Op_SlackAbs) instead of computing LP2(|slew|) itself, which makes that term linear.
 *
 * Only term 2: the raw slew's sign flips across CG steps, the filtered slew's (terms 1/3) does not, and
 * Op_SAFE::signed_terms13 handles those. A separate class so Op_SAFE is untouched; Ax keeps the same layout,
 * so prox() and the root-sum-square limit are inherited.
 */

#include <string>

#include "op_safe.hpp"

namespace Gropt {

// The shared u >= |slew| block: the slew bound is the waveform's, whichever model reads it.
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

    virtual void declare_aux(ProblemData &pd) override;
    virtual void init() override;

    virtual void forward(Eigen::VectorXd &X, Eigen::VectorXd &out) override;
    virtual void transpose(Eigen::VectorXd &X, Eigen::VectorXd &out) override;
    virtual void check(Eigen::VectorXd &X) override;
};

} // namespace Gropt

#endif
