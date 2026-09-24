#include <stdexcept>

#include "spdlog/spdlog.h"

#include "op_slack_abs.hpp"

namespace Gropt {

Op_SlackAbs::Op_SlackAbs(const ProblemData &_pdata, const std::string &_aux_name, double _weight_mod)
    : Operator(_pdata) {
    name = "SlackAbs";
    uses_aux = true;
    aux_name = _aux_name;
    weight_mod = _weight_mod;
}

void Op_SlackAbs::declare_aux(ProblemData &pd) { pd.add_aux(aux_name, pd.n_wave()); }

void Op_SlackAbs::init_aux(Eigen::VectorXd &X_full) {
    // Start tight: u = |dx| satisfies both rows with equality, which is where the solution sits.
    const int off = pdata->aux_offset(aux_name);
    if (off < 0) return;
    for (int j = 0; j < pdata->Naxis; j++) {
        for (int i = 0; i < pdata->N; i++) {
            const int k = j * pdata->N + i;
            const double d = (i == 0) ? X_full(k) : (X_full(k) - X_full(k - 1));
            X_full(off + k) = std::abs(d);
        }
    }
}

void Op_SlackAbs::init() {
    u_offset = pdata->aux_offset(aux_name);
    if (u_offset < 0) {
        throw std::runtime_error("Op_SlackAbs::init: auxiliary block '" + aux_name + "' was never declared");
    }

    target = 0.0;
    tol0 = 0.0;
    tol = 0.0;

    Ax_size = 2 * pdata->Naxis * pdata->N;

    Operator::init();

    // Linear but not trivially normalized: D contributes ~2/dt. Measure it like Op_SAFE does.
    spec_norm = estimate_self_spec_norm(30);
    spec_norm2 = spec_norm * spec_norm;

    spdlog::trace("Op_SlackAbs::init  spec_norm = {:.4f}", spec_norm);
}

void Op_SlackAbs::forward(Eigen::VectorXd &X, Eigen::VectorXd &out) {
    // out = [u - dx ; u + dx], per axis, all in waveform units
    for (int j = 0; j < Naxis; j++) {
        for (int i = 0; i < N; i++) {
            const int k = j * N + i;
            const double d = (i == 0) ? X(k) : (X(k) - X(k - 1));
            const double u = X(u_offset + k);
            out(2 * j * N + i) = u - d;
            out(2 * j * N + i + N) = u + d;
        }
    }
}

void Op_SlackAbs::transpose(Eigen::VectorXd &X, Eigen::VectorXd &out) {
    // Adjoint of [u - dx ; u + dx]: u gets p + q, x gets the difference adjoint of (q - p).
    out.setZero();
    for (int j = 0; j < Naxis; j++) {
        // slack part
        for (int i = 0; i < N; i++) {
            out(u_offset + j * N + i) = X(2 * j * N + i) + X(2 * j * N + i + N);
        }
        // waveform part: form (q - p), then apply D^T, which matches Op_SAFE's differentiator adjoint
        for (int i = 0; i < N; i++) {
            out(j * N + i) = X(2 * j * N + i + N) - X(2 * j * N + i);
        }
        for (int i = 0; i < N - 1; i++) {
            out(j * N + i) = out(j * N + i) - out(j * N + i + 1);
        }
    }
}

void Op_SlackAbs::prox(Eigen::VectorXd &X) {
    if (do_equil) {
        X.array() /= eq_rows.array();
    }
    X.array() *= spec_norm;

    // Both rows are one-sided: project onto the non-negative orthant.
    X = X.cwiseMax(0.0);

    if (do_equil) {
        X.array() *= eq_rows.array();
    }
    X.array() /= spec_norm;
}

void Op_SlackAbs::check(Eigen::VectorXd &X) {
    if (do_equil) {
        X.array() /= eq_rows.array();
    }
    X.array() *= spec_norm;

    // Relative tolerance: the rows carry slew-sized numbers, so an absolute epsilon would be meaningless.
    const double scale = X.cwiseAbs().maxCoeff();
    const int is_feas = (X.minCoeff() >= -1e-8 * std::max(scale, 1.0)) ? 1 : 0;

    if (do_equil) {
        X.array() *= eq_rows.array();
    }
    X.array() /= spec_norm;

    hist_feas.push_back(is_feas);
}

} // namespace Gropt
