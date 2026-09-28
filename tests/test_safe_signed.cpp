// Op_SAFE::signed_terms13: terms 1 and 3 emitted signed, their |.| taken in prox() instead. Their abs is
// outside the low-pass, so the projection takes it exactly (term 2's is inside the filter, so it cannot).
// The set is symmetric in every coordinate, so the projection keeps each sign and only shrinks magnitudes.

#include "op_safe.hpp"
#include "test_util.hpp"

#include <random>

using namespace Gropt;
using namespace gropt_test;

int run_safe_signed_tests() {
    std::printf("\nOp_SAFE signed terms 1 and 3\n");
    int failures = 0;

    const int N = 64, Naxis = 3;
    ProblemData p = make_pdata(N, Naxis, 20e-6, /*pin_ends=*/true);

    Op_SAFE plain(p, 1.0, 1.0);
    plain.safe_params.set_demo_params();
    plain.init();
    Op_SAFE sgn(p, 1.0, 1.0);
    sgn.safe_params.set_demo_params();
    sgn.signed_terms13 = true;
    sgn.init();

    const double sn = plain.spec_norm;               // prox works in the normalized space
    const double upper = 1.0 - plain.cushion;
    const int nt = plain.n_terms;
    std::mt19937 rng(31u);
    std::uniform_real_distribution<double> pos(0.5, 1.5);

    // (1) POSITIVE input, MILD overshoot so no term is driven to the floor: the magnitude-space
    // projection must reproduce the inherited one exactly -- same set, same equal-share step.
    {
        Eigen::VectorXd phys(plain.Ax_size);
        std::uniform_real_distribution<double> tight(0.95, 1.05);
        for (int k = 0; k < phys.size(); k++) phys(k) = tight(rng);
        phys *= 1.02 * upper / (sqrt((double)Naxis) * nt); // just outside the combined limit
        Eigen::VectorXd za = phys / sn, zb = phys / sn;
        plain.prox(za);
        sgn.prox(zb);
        failures += report(za.minCoeff() > 0.0, "mild overshoot clamps nothing", za.minCoeff() * sn);
        failures += report((za - zb).cwiseAbs().maxCoeff() <= 1e-12 * std::max(1.0, za.cwiseAbs().maxCoeff()),
                           "signed prox == inherited prox when nothing clamps",
                           (za - zb).cwiseAbs().maxCoeff());
    }

    // (2) HEAVY overshoot: the inherited prox subtracts a flat share and drives small terms NEGATIVE,
    // outside the non-negative set it is meant to project onto. The magnitude version floors them.
    {
        Eigen::VectorXd phys(plain.Ax_size);
        for (int k = 0; k < phys.size(); k++) phys(k) = pos(rng);
        phys *= upper / phys.maxCoeff();
        Eigen::VectorXd za = phys / sn, zb = phys / sn;
        plain.prox(za);
        sgn.prox(zb);
        failures += report(za.minCoeff() < 0.0, "heavy overshoot: inherited prox goes negative",
                           za.minCoeff() * sn);
        failures += report(zb.minCoeff() >= -1e-15, "heavy overshoot: magnitude prox floors at zero",
                           zb.minCoeff() * sn);
    }

    // (3) SIGNED input: land on the boundary measured with |.|, keep every sign, and shrink every
    // unclamped term of an axis by the same amount (the KKT equal share).
    {
        Eigen::VectorXd vs(plain.Ax_size);
        for (int k = 0; k < vs.size(); k++) vs(k) = pos(rng) * ((k % 3 == 0) ? -1.0 : 1.0);
        vs *= 1.5;
        Eigen::VectorXd zs = vs / sn;
        sgn.prox(zs);
        zs *= sn;

        double worst_rss = 0.0, worst_sign = 0.0, worst_share = 0.0;
        for (int i = 0; i < N; i++) {
            double ss = 0.0;
            for (int j = 0; j < Naxis; j++) {
                double s = 0.0, share = -1.0;
                for (int k = 0; k < nt; k++) {
                    const int idx = j * nt * N + i + k * N;
                    const double zz = zs(idx), vv = vs(idx);
                    if (zz != 0.0 && vv != 0.0 && (zz > 0) != (vv > 0)) worst_sign = 1.0;
                    s += std::abs(zz);
                    if (std::abs(zz) > 1e-9) {
                        const double d = std::abs(vv) - std::abs(zz);
                        if (share < 0.0) share = d;
                        else worst_share = std::max(worst_share, std::abs(d - share));
                    }
                }
                ss += s * s; // stim_thresh_vec is 1 here
            }
            worst_rss = std::max(worst_rss, std::sqrt(ss) - upper);
        }
        failures += report(worst_rss <= 1e-9, "signed prox lands on the limit", worst_rss);
        failures += report(worst_sign == 0.0, "signed prox preserves every sign", worst_sign);
        failures += report(worst_share <= 1e-9, "equal share across an axis's unclamped terms", worst_share);
    }

    // (4) A point already inside the limit must be left alone.
    {
        Eigen::VectorXd phys(plain.Ax_size);
        for (int k = 0; k < phys.size(); k++) phys(k) = pos(rng) * ((k % 2) ? -1.0 : 1.0);
        phys *= 0.1 * upper / (sqrt((double)Naxis) * nt);
        Eigen::VectorXd z = phys / sn;
        sgn.prox(z);
        failures += report((z - phys / sn).cwiseAbs().maxCoeff() == 0.0,
                           "a feasible point is untouched", (z - phys / sn).cwiseAbs().maxCoeff());
    }

    // (5) Scoring must see |.| too: SAFE is even in G, so g and -g score alike, and signed scores like plain.
    {
        Eigen::VectorXd g = Eigen::VectorXd::Zero(N * Naxis);
        for (int j = 0; j < Naxis; j++) {
            for (int i = 1; i < N - 1; i++) {
                g(j * N + i) = 0.04 * std::sin(3.0 * PI * i / (N - 1)) * (1.0 + 0.3 * j);
            }
        }
        for (int it = 0; it < 40 && plain.constraint_violation(g) < 0.5; it++) g *= 2.0; // well over the limit
        const Eigen::VectorXd gn = -g;
        const double vp = plain.constraint_violation(g), vpn = plain.constraint_violation(gn);
        const double vs = sgn.constraint_violation(g), vsn = sgn.constraint_violation(gn);
        failures += report(vp > 0.1, "test waveform violates the limit", vp);
        failures += report(std::abs(vp - vpn) <= 1e-12 * vp, "plain: violation(g) == violation(-g)", vp - vpn);
        failures += report(std::abs(vs - vp) <= 1e-12 * vp, "signed: violation(g) == plain", vs - vp);
        failures += report(std::abs(vsn - vp) <= 1e-12 * vp, "signed: violation(-g) == plain", vsn - vp);

        auto feas = [](Op_SAFE &op, const Eigen::VectorXd &x) {
            Eigen::VectorXd xx = x, ax(op.Ax_size);
            op.forward_op(xx, ax);
            op.check(ax);
            return op.hist_feas.back();
        };
        const int feas_p = feas(plain, g), feas_s = feas(sgn, g), feas_sn = feas(sgn, gn);
        failures += report(feas_p == 0 && feas_s == 0 && feas_sn == 0, "check(): violating g and -g are infeasible",
                           feas_p + feas_s + feas_sn);
    }

    return failures;
}
