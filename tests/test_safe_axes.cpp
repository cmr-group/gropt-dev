// Op_SAFE multi-axis combination: per sample, sqrt(sum_j (stim_j / thresh_j)^2) <= 1, as the scanner adds
// simultaneous axes. Naxis == 1 must reduce to the single-axis clamp exactly.

#include "op_safe.hpp"
#include "test_util.hpp"

#include <memory>

using namespace Gropt;
using namespace gropt_test;

namespace {

// A waveform with enough slew to push SAFE past the threshold on every axis.
Eigen::VectorXd make_wave(int N, int Naxis, bool identical_axes) {
    Eigen::VectorXd x(N * Naxis);
    for (int j = 0; j < Naxis; j++) {
        double amp = identical_axes ? 0.03 : 0.03 * (1.0 + 0.4 * j);
        double frq = identical_axes ? 9.0 : 9.0 + 3.0 * j;
        for (int i = 0; i < N; i++) x(j * N + i) = amp * std::sin(frq * 2.0 * PI * i / N);
    }
    return x;
}

// Peak per-axis stimulation and peak root-sum-square, from a (normalized) Ax.
void stim_stats(Op_SAFE &op, const Eigen::VectorXd &Ax_norm, double &max_axis, double &max_rss) {
    const int N = op.N, Naxis = op.Naxis, nt = op.n_terms;
    Eigen::VectorXd Ax = Ax_norm * op.spec_norm;
    max_axis = 0.0;
    max_rss = 0.0;
    for (int i = 0; i < N; i++) {
        double ss = 0.0;
        for (int j = 0; j < Naxis; j++) {
            double v = Ax(j * nt * N + i) + Ax(j * nt * N + i + N);
            if (nt == 3) v += Ax(j * nt * N + i + 2 * N);
            max_axis = std::max(max_axis, std::abs(v));
            ss += v * v;
        }
        max_rss = std::max(max_rss, std::sqrt(ss));
    }
}

// same_params: give every axis the x-axis coefficients.
std::unique_ptr<Op_SAFE> make_op(ProblemData &p, bool same_params = false) {
    auto op = std::make_unique<Op_SAFE>(p, /*stim_thresh=*/1.0, /*weight_mod=*/1.0);
    op->safe_params.set_demo_params();
    if (same_params) {
        SAFEParams &sp = op->safe_params;
        for (int j = 1; j < 3; j++) {
            sp.tau1[j] = sp.tau1[0]; sp.tau2[j] = sp.tau2[0]; sp.tau3[j] = sp.tau3[0];
            sp.a1[j] = sp.a1[0];     sp.a2[j] = sp.a2[0];     sp.a3[j] = sp.a3[0];
            sp.stim_limit[j] = sp.stim_limit[0];
            sp.g_scale[j] = sp.g_scale[0];
        }
    }
    op->init();
    return op;
}

// Rotate a (Naxis=3) waveform about the z axis then the x axis.
Eigen::VectorXd rotate(const Eigen::VectorXd &x, int N, double az, double ax) {
    Eigen::Matrix3d Rz, Rx;
    Rz << std::cos(az), -std::sin(az), 0, std::sin(az), std::cos(az), 0, 0, 0, 1;
    Rx << 1, 0, 0, 0, std::cos(ax), -std::sin(ax), 0, std::sin(ax), std::cos(ax);
    Eigen::Matrix3d R = Rx * Rz;
    Eigen::VectorXd y(x.size());
    for (int i = 0; i < N; i++) {
        Eigen::Vector3d v(x(i), x(N + i), x(2 * N + i));
        Eigen::Vector3d w = R * v;
        y(i) = w(0); y(N + i) = w(1); y(2 * N + i) = w(2);
    }
    return y;
}

} // namespace

int run_safe_axes_tests() {
    int failures = 0;
    std::printf("\n=== Op_SAFE multi-axis combination ===\n");
    const int N = 400;
    ProblemData p3 = make_pdata(N, 3);

    // 1. Naxis == 1 reduces to the single-axis clamp: prox lands exactly on (1 - cushion) * threshold.
    {
        ProblemData p1 = make_pdata(N, 1);
        Eigen::VectorXd x = make_wave(N, 1, false);
        auto op = make_op(p1);
        Eigen::VectorXd a(op->Ax_size);
        op->forward_op(x, a);
        double ax0, rss0;
        stim_stats(*op, a, ax0, rss0);
        op->prox(a);
        double ax, rss;
        stim_stats(*op, a, ax, rss);
        failures += report(ax0 > 1.0, "Naxis=1: test waveform is over the limit before prox", ax0);
        failures += report(std::abs(ax - (1.0 - op->cushion)) < 1e-9,
                           "Naxis=1: prox lands on (1 - cushion) * threshold", ax);
    }

    // 2. THE REGRESSION: a waveform whose every axis sits exactly at its own limit is NOT feasible, because
    // the scanner adds the axes. Scaling works because SAFE is positively homogeneous of degree 1 in g.
    {
        auto op = make_op(p3);
        Eigen::VectorXd x = make_wave(N, 3, false);
        Eigen::VectorXd a(op->Ax_size);
        op->forward_op(x, a);
        double ax, rss;
        stim_stats(*op, a, ax, rss);
        Eigen::VectorXd xs = x / ax; // every axis now peaks at exactly 1.0
        op->forward_op(xs, a);
        stim_stats(*op, a, ax, rss);
        failures += report(std::abs(ax - 1.0) < 1e-9, "per-axis exactly at the limit", ax);
        failures += report(rss > 1.0, "...yet the combined rss exceeds it (this is the bug being fixed)", rss);
        op->check(a);
        failures += report(op->hist_feas.back() == 0, "check() rejects it", rss);
        failures += report(op->constraint_violation(xs) > 0.0, "constraint_violation reports the excess",
                           op->constraint_violation(xs));
    }

    // 3. prox projects onto the combined limit and its output passes check.
    {
        auto op = make_op(p3);
        Eigen::VectorXd x = make_wave(N, 3, false);
        Eigen::VectorXd a(op->Ax_size);
        op->forward_op(x, a);
        op->prox(a);
        double ax, rss;
        stim_stats(*op, a, ax, rss);
        failures += report(rss <= 1.0, "max rss after prox <= 1", rss);
        failures += report(std::abs(rss - (1.0 - op->cushion)) < 1e-6, "...and lands on the cushion", rss);
        op->check(a);
        failures += report(op->hist_feas.back() == 1, "prox output passes check", rss);
    }

    // 4. Three identical axes: the rss is sqrt(3) x one axis, so each axis is held to 1/sqrt(3).
    {
        Eigen::VectorXd x = make_wave(N, 3, true);
        auto op = make_op(p3, /*same_params=*/true); // equal coefficients too, or the axes are not identical
        Eigen::VectorXd a(op->Ax_size);
        op->forward_op(x, a);
        op->prox(a);
        double ax, rss;
        stim_stats(*op, a, ax, rss);
        failures += report(std::abs(ax - 1.0 / std::sqrt(3.0)) < 0.02, "identical axes: each held to 1/sqrt(3)",
                           ax);
    }

    // 5. SAFE is rotation variant: each axis sums its |.| terms before the axes combine, so equal per-axis
    // coefficients shrink the dependence but cannot remove it.
    {
        Eigen::VectorXd x = make_wave(N, 3, false);
        Eigen::VectorXd xr = rotate(x, N, 0.7, 0.4);
        double rel[2] = {0.0, 0.0};
        for (int same = 0; same < 2; same++) {
            auto op = make_op(p3, same == 1);
            Eigen::VectorXd a(op->Ax_size), b(op->Ax_size);
            op->forward_op(x, a);
            op->forward_op(xr, b);
            double axa, rssa, axb, rssb;
            stim_stats(*op, a, axa, rssa);
            stim_stats(*op, b, axb, rssb);
            rel[same] = std::abs(rssa - rssb) / rssa;
        }
        failures += report(rel[0] > 1e-3, "per-axis coefficients: rotation changes the rss (variant)", rel[0]);
        failures += report(rel[1] > 1e-9, "equal coefficients: still not rotation invariant", rel[1]);
        failures += report(rel[1] < rel[0], "equal coefficients shrink the rotation dependence", rel[1] / rel[0]);
    }
    return failures;
}
