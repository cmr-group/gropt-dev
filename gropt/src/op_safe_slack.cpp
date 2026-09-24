#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "spdlog/spdlog.h"

#include "op_safe_slack.hpp"

namespace Gropt {

Op_SAFE_Slack::Op_SAFE_Slack(const ProblemData &_pdata, double _stim_thresh, double _weight_mod,
                             const std::string &_aux_name)
    : Op_SAFE(_pdata, _stim_thresh, _weight_mod) {
    name = "SAFE_Slack";
    uses_aux = true;
    aux_name = _aux_name;
}

Op_SAFE_Slack::Op_SAFE_Slack(const ProblemData &_pdata, const Eigen::VectorXd &_stim_thresh_vec,
                             double _weight_mod, const std::string &_aux_name)
    : Op_SAFE(_pdata, _stim_thresh_vec, _weight_mod) {
    name = "SAFE_Slack";
    uses_aux = true;
    aux_name = _aux_name;
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

    // term 2 = LP2(u / dt): no abs and no sign. That is the whole point of the lifting -- it makes the
    // operator exactly linear in u, so this block is never re-linearized. The 1/dt converts the slack
    // from waveform units (see Op_SlackAbs) back to the slew the SAFE model expects.
    stim2 = X.segment(u_offset, Ntot) / dt;
    lowpass(stim2, safe_params.alpha2);

    if (n_terms == 3) {
        stim3 = slew_temp;
        lowpass(stim3, safe_params.alpha3);
        if (!signed_terms13) take_abs(stim3, signs3);
    }

    // Same Ax layout as Op_SAFE, so the inherited prox and combined limit apply unchanged.
    for (int j = 0; j < Naxis; j++) {
        for (int i = 0; i < N; i++) {
            out(j * n_terms * N + i) =
                safe_params.a1[j] * stim1(j * N + i) / safe_params.stim_limit[j] * safe_params.g_scale[j];
            out(j * n_terms * N + i + N) =
                safe_params.a2[j] * stim2(j * N + i) / safe_params.stim_limit[j] * safe_params.g_scale[j];
            if (n_terms == 3) {
                out(j * n_terms * N + i + 2 * N) =
                    safe_params.a3[j] * stim3(j * N + i) / safe_params.stim_limit[j] * safe_params.g_scale[j];
            }
        }
    }
}

void Op_SAFE_Slack::transpose(Eigen::VectorXd &X, Eigen::VectorXd &out) {
    out.setZero();

    for (int j = 0; j < Naxis; j++) {
        for (int i = 0; i < N; i++) {
            stim1(j * N + i) =
                safe_params.a1[j] * X(j * n_terms * N + i) / safe_params.stim_limit[j] * safe_params.g_scale[j];
            stim2(j * N + i) = safe_params.a2[j] * X(j * n_terms * N + i + N) / safe_params.stim_limit[j] *
                               safe_params.g_scale[j];
            if (n_terms == 3) {
                stim3(j * N + i) = safe_params.a3[j] * X(j * n_terms * N + i + 2 * N) /
                                   safe_params.stim_limit[j] * safe_params.g_scale[j];
            }
        }
    }

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

void Op_SAFE_Slack::true_Ax(Eigen::VectorXd &out) {
    // Op_SAFE::forward is the unlifted model: its term 2 is LP2(abs(slew)). Run it on the cached waveform
    // with the abs live, and put the captured signs back so the CG's frozen linearization is undisturbed.
    Eigen::VectorXd x = Eigen::VectorXd::Zero(pdata->n_total());
    x.head(Ntot) = x_cache;
    const bool saved = freeze_signs;
    const Eigen::VectorXd s1 = signs1, s2 = signs2, s3 = signs3;
    freeze_signs = false;
    Op_SAFE::forward(x, out);
    signs1 = s1;
    signs2 = s2;
    signs3 = s3;
    freeze_signs = saved;
}

void Op_SAFE_Slack::check(Eigen::VectorXd &X) {
    // X is the LIFTED Ax, whose term 2 is LP2(u). u only reaches abs(slew) as the coupling rows converge,
    // so scoring feasibility from it could call a waveform safe while the true model is over the limit.
    // Score the unlifted model instead: what the scanner would compute.
    (void)X;
    Eigen::VectorXd ax(Ax_size);
    true_Ax(ax);

    int is_feas = 1;
    for (int i = 0; i < N; i++) {
        double ss = 0.0;
        for (int j = 0; j < Naxis; j++) {
            double u = axis_stim(ax, j, i) / stim_thresh_vec(j * N + i);
            ss += u * u;
        }
        if (sqrt(ss) > 1.0) is_feas = 0;
    }
    hist_feas.push_back(is_feas);
}

double Op_SAFE_Slack::constraint_violation(const Eigen::VectorXd &x_new) {
    // Same reasoning as check(): report the true model's fractional excess, not the lifted one's.
    const Eigen::VectorXd keep = x_cache;
    x_cache = x_new.head(Ntot);
    Eigen::VectorXd ax(Ax_size);
    true_Ax(ax);
    x_cache = keep;

    double worst = 0.0;
    for (int i = 0; i < N; i++) {
        double ss = 0.0;
        for (int j = 0; j < Naxis; j++) {
            double u = axis_stim(ax, j, i) / stim_thresh_vec(j * N + i);
            ss += u * u;
        }
        worst = std::max(worst, sqrt(ss) - 1.0);
    }
    return worst;
}

} // namespace Gropt
