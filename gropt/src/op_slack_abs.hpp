#ifndef OP_SLACK_ABS_H
#define OP_SLACK_ABS_H

/**
 * u >= |slew|, as the two linear rows u - Dx >= 0 and u + Dx >= 0. Op_SAFE_Slack reads u in place of
 * |slew| so its term 2 is linear; this operator is what keeps u honest. Tight without extra machinery:
 * every coefficient that sees u is non-negative, so the only pressure on u is downward.
 *
 * u is held in WAVEFORM units, u >= |x_i - x_{i-1}|, i.e. dt times the slew. Same bound, but the primal
 * is one vector whose proximal term weights every entry alike, and slew-sized numbers (~200) beside a
 * waveform of ~0.08 are 2500x out of scale -- the solve drives u to zero, the coupling can then only be
 * met by driving the waveform to zero too, and the reweighting chases it there. Op_SAFE_Slack divides
 * by dt when it reads the block.
 *
 * One instance serves every Op_SAFE_Slack sharing the block, so GroptParams adds it once.
 */

#include "Eigen/Dense"
#include <string>

#include "op_main.hpp"

namespace Gropt {

class Op_SlackAbs : public Operator {
  protected:
    std::string aux_name;
    int u_offset = -1;

  public:
    Op_SlackAbs(const ProblemData &_pdata, const std::string &_aux_name, double _weight_mod);

    const std::string &aux_block() const { return aux_name; }

    virtual void declare_aux(ProblemData &pd) override;
    virtual void init_aux(Eigen::VectorXd &X_full) override;
    virtual void init();

    virtual void forward(Eigen::VectorXd &X, Eigen::VectorXd &out);
    virtual void transpose(Eigen::VectorXd &X, Eigen::VectorXd &out);
    virtual void prox(Eigen::VectorXd &X);
    virtual void check(Eigen::VectorXd &X);

    // Ax is [u - Dx ; u + Dx]: 2 blocks of N per axis.
    virtual std::vector<int> Ax_block_lengths() const override { return std::vector<int>(2 * Naxis, N); }
};

} // namespace Gropt

#endif
