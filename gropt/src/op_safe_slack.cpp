#include <stdexcept>

#include "op_safe_slack.hpp"

namespace Gropt {

Op_SAFE_Slack::Op_SAFE_Slack(const ProblemData &_pdata, double _stim_thresh, double _weight_mod,
                             const std::string &_aux_name)
    : Op_SAFE(_pdata, _stim_thresh, _weight_mod), aux_name(_aux_name) {
    name = "SAFE_Slack";
    uses_aux = true;
}

Op_SAFE_Slack::Op_SAFE_Slack(const ProblemData &_pdata, const Eigen::VectorXd &_stim_thresh_vec,
                             double _weight_mod, const std::string &_aux_name)
    : Op_SAFE(_pdata, _stim_thresh_vec, _weight_mod), aux_name(_aux_name) {
    name = "SAFE_Slack";
    uses_aux = true;
}

void Op_SAFE_Slack::declare_aux(ProblemData &pd) { pd.add_aux(aux_name, pd.n_wave()); }

void Op_SAFE_Slack::init() {
    // Resolve before Op_SAFE::init(), which power-iterates this operator's own forward/transpose.
    u_offset = pdata->aux_offset(aux_name);
    if (u_offset < 0) {
        throw std::runtime_error("Op_SAFE_Slack::init: auxiliary block was never declared: " + aux_name);
    }
    x_cache.setZero(pdata->n_wave());
    Op_SAFE::init();
}

void Op_SAFE_Slack::forward(Eigen::VectorXd &X, Eigen::VectorXd &out) {
    x_cache = X.head(Ntot);
    compute_slew(X);

    stim1 = slew_temp;                                   // terms 1 and 3 exactly as Op_SAFE does them
    lowpass(stim1, safe_params.alpha1);
    if (!signed_terms13) take_abs(stim1, signs1);

    // term 2 = LP2(u / dt): linear in u, never re-linearized; 1/dt turns the waveform-unit slack into slew
    stim2 = X.segment(u_offset, Ntot) / dt;
    lowpass(stim2, safe_params.alpha2);

    if (n_terms == 3) {
        stim3 = slew_temp;
        lowpass(stim3, safe_params.alpha3);
        if (!signed_terms13) take_abs(stim3, signs3);
    }

    emit(out); // same Ax layout as Op_SAFE
}

void Op_SAFE_Slack::transpose(Eigen::VectorXd &X, Eigen::VectorXd &out) {
    out.setZero();
    absorb(X);

    if (!signed_terms13) stim1.array() *= signs1.array();
    lowpass_T(stim1, safe_params.alpha1);

    // term 2's adjoint lands on the SLACK block -- no sign and no D^T, because u is its input
    lowpass_T(stim2, safe_params.alpha2);
    out.segment(u_offset, Ntot) = stim2 / dt;

    if (n_terms == 3) {
        if (!signed_terms13) stim3.array() *= signs3.array();
        lowpass_T(stim3, safe_params.alpha3);
    }

    // only the filtered-then-abs terms reach the waveform, through D^T
    out.head(Ntot) = (n_terms == 3) ? (stim1 + stim3) : stim1;
    diff_T(out);
}

void Op_SAFE_Slack::check(Eigen::VectorXd &) {
    // Score the unlifted model (Op_SAFE::constraint_violation): LP2(u) only reaches LP2(|slew|) as the
    // coupling converges, so the lifted Ax could call an over-limit waveform safe.
    hist_feas.push_back(constraint_violation(x_cache) > 0.0 ? 0 : 1);
}

} // namespace Gropt
