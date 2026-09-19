#ifndef STEP_MONITOR_H
#define STEP_MONITOR_H

/**
 * Trust-region step acceptance for the SDMM solver: on reject, the solver re-solves from the previous
 * iterate with a larger proximal sigma.
 */

#include <memory>
#include <string>

#include "Eigen/Dense"

namespace Gropt {

class GroptParams;

struct StepDecision {
    bool accept = true;       // accept the CG step, or reject and re-solve with a larger sigma
    double sigma_scale = 1.0; // multiply the proximal sigma by this on reject (>1); ignored on accept
    double signal = 0.0;      // monitored value, for logging
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

// Default monitor: rejects when the largest Operator::linearization_error(x_new) exceeds tol, i.e. the step
// left a nonlinear op's frozen-linearization region (e.g. SAFE's |.|). Linear ops report 0.
class LinearizationErrorMonitor : public StepMonitor {
  public:
    StepDecision check(GroptParams &gp, const Eigen::VectorXd &x_old, const Eigen::VectorXd &x_new,
                       double cur_sigma) override;
    const char *name() const override { return "linearization_error"; }
};

// Rejects when the relative step ||x_new - x_old|| / ||x_old|| exceeds tol. Operator-agnostic, but the tol
// needs tuning.
class RelStepMonitor : public StepMonitor {
  public:
    RelStepMonitor() { tol = 0.5; }
    StepDecision check(GroptParams &gp, const Eigen::VectorXd &x_old, const Eigen::VectorXd &x_new,
                       double cur_sigma) override;
    const char *name() const override { return "rel_step"; }
};

// Feasibility funnel: rejects when the constraint violation (sum over ops of
// Operator::constraint_violation) grows beyond max(previous, tol).
class FeasibilityMonitor : public StepMonitor {
  public:
    FeasibilityMonitor() { tol = 0.02; }
    StepDecision check(GroptParams &gp, const Eigen::VectorXd &x_old, const Eigen::VectorXd &x_new,
                       double cur_sigma) override;
    const char *name() const override { return "feasibility"; }
};

// Monitor by name ("linearization_error" | "rel_step" | "feasibility"); nullptr for "none"/unknown.
// tol/bump override the defaults only when > 0.
std::unique_ptr<StepMonitor> make_step_monitor(const std::string &name, double tol, double bump);

} // namespace Gropt
  
#endif
