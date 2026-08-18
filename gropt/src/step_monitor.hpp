#ifndef STEP_MONITOR_H
#define STEP_MONITOR_H

/**
 * Pluggable trust-region step-acceptance strategies for the SDMM solver. After the inner CG produces a
 * candidate iterate, a StepMonitor decides whether the linearized model held over the step (accept) or
 * was outrun (reject -> the solver re-solves from the previous iterate with the proximal sigma scaled
 * up).
 */

#include <memory>
#include <string>

#include "Eigen/Dense"

namespace Gropt {

class GroptParams; // forward decl (check() reads gp.all_op)

// Outcome of a step-acceptance test.
struct StepDecision {
    bool accept = true;       // accept the CG step, or reject and re-solve with a larger sigma
    double sigma_scale = 1.0; // multiply the proximal sigma by this on reject (>1); ignored on accept
    double signal = 0.0;      // the raw signal value (for logging / diagnostics)
};

// Base strategy. tol = reject threshold (meaning depends on the signal); bump = sigma multiplier on reject.
class StepMonitor {
  public:
    virtual ~StepMonitor() = default;
    virtual StepDecision check(GroptParams &gp, const Eigen::VectorXd &x_old,
                               const Eigen::VectorXd &x_new, double cur_sigma) = 0;
    virtual const char *name() const = 0;

    double tol = 0.2;
    double bump = 4.0;
};

// DEFAULT: model fidelity of every nonlinear operator's frozen linearization (gain-ratio).
// Rejects when max over ops of Operator::linearization_error(x_new) exceeds tol -- for SAFE this is
// exactly "did the step leave the |.| linearization's valid region"; linear ops report 0 (never trigger).
class LinearizationErrorMonitor : public StepMonitor {
  public:
    StepDecision check(GroptParams &gp, const Eigen::VectorXd &x_old, const Eigen::VectorXd &x_new,
                       double cur_sigma) override;
    const char *name() const override { return "linearization_error"; }
};

// ALTERNATIVE: relative primal step ||x_new - x_old|| / ||x_old|| (the "change in x" divergence signal).
// Operator-agnostic and cheap, but needs a tuned tol (no self-calibration like the model-fidelity one).
class RelStepMonitor : public StepMonitor {
  public:
    RelStepMonitor() { tol = 0.5; }
    StepDecision check(GroptParams &gp, const Eigen::VectorXd &x_old, const Eigen::VectorXd &x_new,
                       double cur_sigma) override;
    const char *name() const override { return "rel_step"; }
};

// FEASIBILITY FUNNEL (SQP restoration): accept a step only if the worst-sample constraint violation
// does not grow beyond max(previous, tol).
class FeasibilityMonitor : public StepMonitor {
  public:
    FeasibilityMonitor() { tol = 0.02; }
    StepDecision check(GroptParams &gp, const Eigen::VectorXd &x_old, const Eigen::VectorXd &x_new,
                       double cur_sigma) override;
    const char *name() const override { return "feasibility"; }
};

// By name ("linearization_error" | "rel_step" | "feasibility"); nullptr for "none"/unknown (TR
// off). tol / bump override the monitor's own default only when > 0 (pass tr_tol<=0 to keep the default,
// since the good tol differs per monitor: ~0.2 linearization, ~0.02 feasibility).
std::unique_ptr<StepMonitor> make_step_monitor(const std::string &name, double tol, double bump);

} // namespace Gropt
  
#endif
