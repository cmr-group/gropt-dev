// Raster conversion: resampling onto a raster that divides the solve raster is exact. Pins that, plus the
// quadratures changed so it holds (Op_Concomitant::exact_quad, SAFEParams::alpha_exact, pwl moments).

#include "gropt_utils.hpp"
#include "op_bvalue.hpp"
#include "op_concomitant.hpp"
#include "raster.hpp"
#include "test_util.hpp"

#include <random>

using namespace Gropt;
using namespace gropt_test;

namespace {

// A spin-echo-shaped single-axis test waveform on dt: zero at both ends and across a 180 window, smooth
// in between, sign-flipped after the 180 -- the layout diff_init() produces.
Eigen::VectorXd make_wave(int N, int i180s, int i180e) {
    Eigen::VectorXd g = Eigen::VectorXd::Zero(N);
    for (int i = 1; i < i180s; i++) {
        const double u = static_cast<double>(i) / i180s;
        g(i) = 0.05 * std::sin(PI * u) * (1.0 + 0.3 * std::cos(5.0 * PI * u));
    }
    for (int i = i180e + 1; i < N - 1; i++) {
        const double u = static_cast<double>(i - i180e) / (N - 1 - i180e);
        g(i) = -0.05 * std::sin(PI * u) * (1.0 + 0.2 * std::cos(3.0 * PI * u));
    }
    return g;
}

Eigen::VectorXd make_inv(int N, int i_inv) {
    Eigen::VectorXd v = Eigen::VectorXd::Ones(N);
    for (int i = i_inv; i < N; i++) v(i) = -1.0;
    return v;
}

int test_raster_ratio() {
    std::printf("\n-- raster_ratio --\n");
    int f = 0;
    f += check_close(raster_ratio(400e-6, 10e-6), 40, 0, "400 us / 10 us == 40");
    f += check_close(raster_ratio(10e-6, 10e-6), 1, 0, "10 us / 10 us == 1");
    for (auto bad : {std::pair<double, double>{25e-6, 10e-6},  // 2.5, not an integer
                     std::pair<double, double>{10e-6, 400e-6}, // refining only
                     std::pair<double, double>{0.0, 10e-6}}) {
        f += report(throws_invalid([&] { raster_ratio(bad.first, bad.second); }),
                    "rejects dt_src=" + std::to_string(bad.first) + " dt_tgt=" + std::to_string(bad.second));
    }
    return f;
}

int test_resample() {
    std::printf("\n-- resample_waveform --\n");
    int f = 0;
    const int N = 117, R = 40;
    const double dt = 400e-6, dt_f = 10e-6;
    const Eigen::VectorXd g = make_wave(N, 71, 86);
    const Eigen::VectorXd gf = resample_waveform(g, 1, dt, dt_f);

    f += check_close(gf.size(), (N - 1) * R + 1, 0, "length == (N-1)*R + 1");

    double node_err = 0.0;
    for (int j = 0; j < N; j++) node_err = std::max(node_err, std::abs(gf(j * R) - g(j)));
    f += check_close(node_err, 0.0, 0.0, "source samples reproduced exactly");

    f += check_close(gf.cwiseAbs().maxCoeff(), g.cwiseAbs().maxCoeff(), 1e-15, "gmax unchanged");

    double s_c = 0.0, s_f = 0.0;
    for (int i = 1; i < N; i++) s_c = std::max(s_c, std::abs(g(i) - g(i - 1)) / dt);
    for (int i = 1; i < gf.size(); i++) s_f = std::max(s_f, std::abs(gf(i) - gf(i - 1)) / dt_f);
    f += check_close(s_f, s_c, 1e-9, "smax unchanged");

    const Eigen::VectorXd padded = resample_waveform(g, 1, dt, dt_f, static_cast<int>(gf.size()) + 25);
    f += report(padded.size() == gf.size() + 25 && padded.tail(25).cwiseAbs().maxCoeff() == 0.0 &&
                    padded.head(gf.size()).isApprox(gf),
                "N_out pads with zeros and leaves the waveform alone");
    f += report(throws_invalid([&] { resample_waveform(g, 1, dt, dt_f, static_cast<int>(gf.size()) - 1); }),
                "N_out below the natural length is rejected");

    // Two axes stay independent (axis 1 is a scaled copy, so resampling must scale the same way)
    Eigen::VectorXd g2(2 * N);
    g2.head(N) = g;
    g2.tail(N) = -2.0 * g;
    const Eigen::VectorXd g2f = resample_waveform(g2, 2, dt, dt_f);
    const int Nf = static_cast<int>(gf.size());
    f += check_close((g2f.head(Nf) - gf).cwiseAbs().maxCoeff(), 0.0, 1e-15, "2-axis: axis 0 matches");
    f += check_close((g2f.tail(Nf) + 2.0 * gf).cwiseAbs().maxCoeff(), 0.0, 1e-15, "2-axis: axis 1 matches");
    return f;
}

// Reference quadratures mirroring Op_Moment / Op_BValue, kept in the test so the library ships only one.

// Op_Moment: 1e6 * dt * (1000 dt i)^k, or the tent weight when pwl
double ref_moment(const Eigen::VectorXd &g, double dt, int k, bool pwl) {
    const double h_ms = 1000.0 * dt;
    double m = 0.0;
    for (int i = 0; i < g.size(); i++) {
        const double T = h_ms * i;
        const double w = pwl ? pwl_moment_weight(k, T, h_ms, i > 0, i < g.size() - 1)
                             : h_ms * std::pow(T, static_cast<double>(k));
        m += 1000.0 * w * g(i);
    }
    return m;
}

// Op_BValue: inclusive running sum, rectangle sum of its square
double ref_bvalue(const Eigen::VectorXd &g, double dt, const Eigen::VectorXd &iv) {
    const double scale = std::pow(267.5221900e6 / 1000.0 * dt, 2) * dt;
    double gt = 0.0, b = 0.0;
    for (int i = 0; i < g.size(); i++) {
        gt += g(i) * iv(i);
        b += gt * gt * scale;
    }
    return b;
}

int test_bvalue_drift() {
    std::printf("\n-- b-value across rasters --\n");
    int f = 0;
    const int N = 117, i180s = 71, i180e = 86, i_inv = 78;
    const double dt = 400e-6, dt_f = 10e-6;
    Eigen::VectorXd g = make_wave(N, i180s, i180e);
    const Eigen::VectorXd iv = make_inv(N, i_inv);

    // Null M0 as a solve would; the b-value drift closed form below needs it.
    double m0 = 0.0, w = 0.0;
    for (int i = 0; i < N; i++) {
        m0 += iv(i) * g(i);
        w += (iv(i) < 0 && g(i) != 0.0) ? 1.0 : 0.0;
    }
    for (int i = 0; i < N; i++) {
        if (iv(i) < 0 && g(i) != 0.0) g(i) += m0 / w; // iv = -1 there, so this subtracts from M0
    }

    const Eigen::VectorXd gf = resample_waveform(g, 1, dt, dt_f);
    const Eigen::VectorXd ivf = resample_inv_vec(iv, 1, raster_ratio(dt, dt_f));

    // b-value drifts by gamma^2 (dt^2 - dt_f^2)/4 * integral(g^2), and by nothing else
    double int_g2 = 0.0;
    for (int i = 0; i + 1 < gf.size(); i++) {
        const double a = gf(i), b = gf(i + 1);
        int_g2 += dt_f / 3.0 * (a * a + a * b + b * b);
    }
    const double gam = 267.5221900e6 / 1000.0;
    const double pred = gam * gam * (dt * dt - dt_f * dt_f) / 4.0 * int_g2;
    const double b_c = ref_bvalue(g, dt, iv), b_f = ref_bvalue(gf, dt_f, ivf);
    f += check_close(b_c - b_f, pred, 2e-2, "b-value drift matches the closed form");
    f += report(std::abs(b_c - b_f) < 0.01 * b_c, "b-value drift under 1%");
    return f;
}

int test_safe_alpha_exact() {
    std::printf("\n-- SAFE alpha_exact raster invariance --\n");
    int f = 0;
    const int N = 117;
    const double dt = 400e-6, dt_f = 10e-6;
    const Eigen::VectorXd g = make_wave(N, 71, 86);
    const Eigen::VectorXd gf = resample_waveform(g, 1, dt, dt_f);

    const double euler_c = get_SAFE_eigen(g, 1, dt, true, 0, false).maxCoeff();
    const double euler_f = get_SAFE_eigen(gf, 1, dt_f, true, 0, false).maxCoeff();
    const double exact_c = get_SAFE_eigen(g, 1, dt, true, 0, true).maxCoeff();
    const double exact_f = get_SAFE_eigen(gf, 1, dt_f, true, 0, true).maxCoeff();
    std::printf("       euler: %.6f (400 us) vs %.6f (10 us);  exact: %.6f vs %.6f\n", euler_c, euler_f, exact_c,
                exact_f);

    f += check_close(exact_f, exact_c, 1e-6, "alpha_exact: same PNS on both rasters");
    f += report(std::abs(euler_f - euler_c) > 10.0 * std::abs(exact_f - exact_c),
                "backward Euler is raster dependent, and much more so than exact");
    f += report(euler_c < exact_c, "backward Euler under-reports at the coarse raster");
    return f;
}

// Op_Concomitant exposes energies() to subclasses; this reaches it for the finite-difference check.
class ConcomitantProbe : public Op_Concomitant {
  public:
    using Op_Concomitant::Op_Concomitant;
    void probe(const Eigen::VectorXd &X, double &pos, double &neg, Eigen::VectorXd *dp, Eigen::VectorXd *dn) const {
        energies(X, pos, neg, dp, dn);
    }
};

int test_concomitant_energies() {
    std::printf("\n-- Op_Concomitant energies and gradient --\n");
    int f = 0;
    const int N = 117, i_inv = 78;
    ProblemData p = make_pdata(N, 1, 400e-6);
    p.inv_vec = make_inv(N, i_inv);

    const Eigen::VectorXd g = make_wave(N, 71, 86);

    for (bool exact : {true, false}) {
        ConcomitantProbe op(p, 0, true, 1.0, 0.1, 1.0);
        op.exact_quad = exact;
        const std::string tag = exact ? "exact_quad" : "legacy";

        double pos = 0.0, neg = 0.0;
        Eigen::VectorXd dp, dn;
        op.probe(g, pos, neg, &dp, &dn);

        // Finite differences on a handful of samples, one per side
        std::mt19937 rng(7);
        std::uniform_int_distribution<int> pre(1, 60), post(95, N - 2);
        double worst = 0.0, scale = 0.0;
        for (int k = 0; k < 8; k++) {
            const int i = (k % 2 == 0) ? pre(rng) : post(rng);
            const double h = 1e-7;
            Eigen::VectorXd gp_ = g, gm = g;
            gp_(i) += h;
            gm(i) -= h;
            double pp, np_, pm, nm;
            op.probe(gp_, pp, np_, nullptr, nullptr);
            op.probe(gm, pm, nm, nullptr, nullptr);
            worst = std::max(worst, std::abs((pp - pm) / (2 * h) - dp(i)));
            worst = std::max(worst, std::abs((np_ - nm) / (2 * h) - dn(i)));
            scale = std::max(scale, std::max(std::abs(dp(i)), std::abs(dn(i))));
        }
        f += report(worst < 1e-6 * std::max(scale, 1e-12), tag + ": gradient matches finite differences");

        // pos and neg must stay functions of disjoint sample sets, or prox's rescaling is invalid
        bool disjoint = true;
        for (int i = 0; i < N; i++) {
            if (dp(i) != 0.0 && dn(i) != 0.0) disjoint = false;
        }
        f += report(disjoint, tag + ": pos and neg depend on disjoint samples");
    }

    // The legacy rectangle sum drifts with the raster; the exact one does not.
    const Eigen::VectorXd gf = resample_waveform(g, 1, 400e-6, 10e-6);
    ProblemData pf = make_pdata(static_cast<int>(gf.size()), 1, 10e-6);
    pf.inv_vec = resample_inv_vec(p.inv_vec, 1, 40);

    for (bool exact : {true, false}) {
        ConcomitantProbe oc(p, 0, true, 1.0, 0.1, 1.0), of(pf, 0, true, 1.0, 0.1, 1.0);
        oc.exact_quad = exact;
        of.exact_quad = exact;
        double pc, nc, pf_, nf;
        oc.probe(g, pc, nc, nullptr, nullptr);
        of.probe(gf, pf_, nf, nullptr, nullptr);
        const double rc = pc / nc, rf = pf_ / nf;
        std::printf("       %-11s ratio %.8f (400 us) vs %.8f (10 us)\n", exact ? "exact_quad" : "legacy", rc, rf);
        if (exact) {
            f += check_close(rf, rc, 1e-9, "exact_quad: ratio invariant across rasters");
        } else {
            f += report(std::abs(rf - rc) > 1e-4 * rc, "legacy: ratio is raster dependent (documented)");
        }
    }
    return f;
}

int test_pwl_moment_weights() {
    std::printf("\n-- pwl_moment_weight --\n");
    int f = 0;
    const double T = 3.7, h = 0.4;

    // Closed forms written out independently of the implementation's loop
    f += check_close(pwl_moment_weight(0, T, h, true, true), h, 1e-12, "k=0 -> h");
    f += check_close(pwl_moment_weight(1, T, h, true, true), h * T, 1e-12, "k=1 -> h T");
    f += check_close(pwl_moment_weight(2, T, h, true, true), h * (T * T + h * h / 6.0), 1e-12,
                     "k=2 -> h (T^2 + h^2/6)");
    f += check_close(pwl_moment_weight(3, T, h, true, true), h * (T * T * T + h * h * T / 2.0), 1e-12,
                     "k=3 -> h (T^3 + h^2 T/2)");
    f += check_close(pwl_moment_weight(4, T, h, true, true),
                     h * (std::pow(T, 4) + h * h * T * T + std::pow(h, 4) / 15.0), 1e-12,
                     "k=4 -> h (T^4 + h^2 T^2 + h^4/15)");

    // Half tents at T=0: right = h^(k+1)/((k+1)(k+2)), left = (-1)^k right (odd k cancel).
    for (int k = 0; k <= 4; k++) {
        const double right = pwl_moment_weight(k, 0.0, h, false, true);
        const double left = pwl_moment_weight(k, 0.0, h, true, false);
        const double full = pwl_moment_weight(k, 0.0, h, true, true);
        const double expect = std::pow(h, k + 1) / ((k + 1.0) * (k + 2.0));
        const double sign = (k % 2 == 0) ? 1.0 : -1.0;
        const std::string K = "k=" + std::to_string(k);
        f += check_close(right, expect, 1e-12, K + ": right half tent at T=0 = h^(k+1)/((k+1)(k+2))");
        f += check_close(left, sign * expect, 1e-12, K + ": left half tent at T=0 = (-1)^k * right");
        f += check_close(full, right + left, 1e-12, K + ": halves sum to the full weight");
    }

    // Integrating t^k against the tent basis must reproduce the integral of a known waveform.
    // g(t) = t on [0, (n-1)h] is exactly representable, so integral t*t^k = L^(k+2)/(k+2).
    const int n = 25;
    for (int k = 0; k <= 4; k++) {
        double m = 0.0;
        for (int i = 0; i < n; i++) {
            m += pwl_moment_weight(k, i * h, h, i > 0, i < n - 1) * (i * h);
        }
        const double L = (n - 1) * h;
        f += check_close(m, std::pow(L, k + 2) / (k + 2), 1e-10,
                         "k=" + std::to_string(k) + ": exact for g(t) = t");
    }
    return f;
}

int test_moment_raster_invariance() {
    std::printf("\n-- moments across rasters, M0 NOT nulled --\n");
    int f = 0;
    const int N = 41, R = 20;
    const double h = 400e-6;
    Eigen::VectorXd g = make_wave(N, 21, 24); // no moment nulling applied, so M0 != 0
    const Eigen::VectorXd gf = resample_waveform(g, 1, R);

    for (bool pwl : {true, false}) {
        for (int k = 0; k <= 3; k++) {
            const double mc = ref_moment(g, h, k, pwl);
            const double mf = ref_moment(gf, h / R, k, pwl);
            const double rel = std::abs(mf - mc) / (std::abs(mc) + 1e-300);
            const bool invariant = pwl || k < 2;
            const std::string nm = "pwl_quad=" + std::to_string(pwl) + " k=" + std::to_string(k) +
                                   (invariant ? " invariant" : " drifts (documented)");
            f += report(invariant ? rel < 1e-10 : rel > 1e-6, nm);
        }
    }
    return f;
}

// Independent closed form for gamma^2 * integral q^2 over the piecewise-linear waveform: on each
// interval q = Q + G u + (S/2) u^2, so the integral is a polynomial in (Q, G, S, h).
double ref_bvalue_pwl(const Eigen::VectorXd &g, double dt, const Eigen::VectorXd &iv) {
    const double gam = 267.5221900e6 / 1000.0;
    double Q = 0.0, tot = 0.0;
    for (int i = 0; i + 1 < g.size(); i++) {
        const double G = g(i) * (iv.size() ? iv(i) : 1.0);
        const double Gn = g(i + 1) * (iv.size() ? iv(i + 1) : 1.0);
        const double S = (Gn - G) / dt;
        tot += Q * Q * dt + Q * G * dt * dt + (G * G / 3.0 + Q * S / 3.0) * dt * dt * dt +
               (G * S / 4.0) * std::pow(dt, 4) + (S * S / 20.0) * std::pow(dt, 5);
        Q += dt * (G + Gn) / 2.0;
    }
    return gam * gam * tot;
}

// Op_BValue's own reported value on a given raster.
double op_bvalue(const Eigen::VectorXd &g, double dt, const Eigen::VectorXd &iv, bool pwl) {
    ProblemData p = make_pdata(static_cast<int>(g.size()), 1, dt);
    p.inv_vec = iv;
    Op_BValue op(p, 1.0, 1.0, -1, -1, 1.0, BVALUE_MODE_MINVAL, 1.01);
    op.pwl_quad = pwl;
    op.init();
    Eigen::VectorXd x = g;
    return op.get_bvalue(x);
}

int test_bvalue_pwl_quad() {
    std::printf("\n-- Op_BValue pwl_quad --\n");
    int f = 0;
    const int N = 117, R = 40;
    const double dt = 400e-6, dt_f = 10e-6;
    Eigen::VectorXd g = make_wave(N, 71, 86);
    const Eigen::VectorXd iv = make_inv(N, 78);
    const Eigen::VectorXd gf = resample_waveform(g, 1, dt, dt_f);
    const Eigen::VectorXd ivf = resample_inv_vec(iv, 1, R);

    // the operator must reproduce the independent closed form, on any raster
    f += check_close(op_bvalue(g, dt, iv, true), ref_bvalue_pwl(g, dt, iv), 1e-10,
                     "pwl value matches the closed form at 400 us");
    f += check_close(op_bvalue(gf, dt_f, ivf, true), ref_bvalue_pwl(gf, dt_f, ivf), 1e-10,
                     "pwl value matches the closed form at 10 us");

    // and it must not depend on the raster, which the default quadrature does
    const double p_c = op_bvalue(g, dt, iv, true), p_f = op_bvalue(gf, dt_f, ivf, true);
    const double r_c = op_bvalue(g, dt, iv, false), r_f = op_bvalue(gf, dt_f, ivf, false);
    std::printf("       pwl  %.6f -> %.6f     default  %.6f -> %.6f\n", p_c, p_f, r_c, r_f);
    f += check_close(p_f, p_c, 1e-9, "pwl b-value invariant across rasters");
    f += report(std::abs(r_f - r_c) > 1e-5 * r_c, "default b-value drifts (documented)");

    // the default sits high of the truth by the closed-form bias, and pwl sits on it
    f += report(r_c > p_c, "default reads high of the exact value at the coarse raster");
    f += report(std::abs(p_c - ref_bvalue_pwl(gf, dt_f, ivf)) < 1e-9 * p_c,
                    "coarse pwl value equals the 10 us truth");
    return f;
}

} // namespace

int run_raster_tests() {
    std::printf("\n=== raster / quadrature tests ===\n");
    int f = 0;
    f += test_raster_ratio();
    f += test_resample();
    f += test_bvalue_drift();
    f += test_safe_alpha_exact();
    f += test_concomitant_energies();
    f += test_pwl_moment_weights();
    f += test_moment_raster_invariance();
    f += test_bvalue_pwl_quad();
    return f;
}
