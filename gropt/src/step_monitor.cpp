#include "step_monitor.hpp"

#include <algorithm>

#include "gropt_params.hpp"
#include "op_main.hpp"

namespace Gropt {

StepDecision LinearizationErrorMonitor::check(GroptParams &gp, const Eigen::VectorXd & /*x_old*/,
                                              const Eigen::VectorXd &x_new, double /*cur_sigma*/) {
    double err = 0.0;
    for (auto &op : gp.all_op) {
        err = std::max(err, op->linearization_error(x_new));
    }
    if (err > tol) {
        return {false, bump, err};
    }
    return {true, 1.0, err};
}

StepDecision RelStepMonitor::check(GroptParams & /*gp*/, const Eigen::VectorXd &x_old,
                                   const Eigen::VectorXd &x_new, double /*cur_sigma*/) {
    double on = x_old.norm();
    double rel = (on > 0.0) ? (x_new - x_old).norm() / on : 0.0;
    if (rel > tol) {
        return {false, bump, rel};
    }
    return {true, 1.0, rel};
}

StepDecision FeasibilityMonitor::check(GroptParams &gp, const Eigen::VectorXd &x_old,
                                       const Eigen::VectorXd &x_new, double /*cur_sigma*/) {
    double v_old = 0.0, v_new = 0.0;
    for (auto &op : gp.all_op) {
        v_old += op->constraint_violation(x_old);
        v_new += op->constraint_violation(x_new);
    }
    // Accept iff the worst-sample violation didn't grow past max(previous, feasibility slack). Rejects the
    // objective creeping past the constraint AND blow-ups (huge violation); the ADMM dual reduces v_old.
    if (v_new <= std::max(v_old, tol)) {
        return {true, 1.0, v_new};
    }
    return {false, bump, v_new};
}

std::unique_ptr<StepMonitor> make_step_monitor(const std::string &name, double tol, double bump) {
    std::unique_ptr<StepMonitor> m;
    if (name == "linearization_error" || name == "linearization" || name == "lin") {
        m = std::make_unique<LinearizationErrorMonitor>();
    } else if (name == "rel_step" || name == "relstep") {
        m = std::make_unique<RelStepMonitor>();
    } else if (name == "feasibility" || name == "feas") {
        m = std::make_unique<FeasibilityMonitor>();
    } else {
        return nullptr; // "none" / unknown -> trust region disabled
    }
    if (tol > 0.0) m->tol = tol; // pass tr_tol<=0 to keep the monitor's own default (differs per monitor)
    if (bump > 0.0) m->bump = bump;
    return m;
}

} // namespace Gropt
