// Op_SAFE::signed_terms13 -- terms 1 and 3 emitted signed, their |.| taken in prox() instead.
//
// Their abs is OUTSIDE the low-pass, so the constraint at a sample depends on only a few per-sample
// coordinates and the projection can take the absolute values exactly. Term 2 cannot be done this way:
// its abs is inside the filter, so LP2(|slew|) at a sample depends on every earlier one.
//
// The set is symmetric in every coordinate, so the projection keeps each sign and only shrinks
// magnitudes. These tests pin that it agrees with the inherited prox wherever the two should agree, and
// that it fixes the one place they differ.

#include "Eigen/Dense"
#include "op_safe.hpp"
#include "problem_data.hpp"

#include <cmath>
#include <cstdio>
#include <random>
#include <string>

using namespace Gropt;

namespace {

int report(bool ok, const std::string &name, double value) {
    std::printf("  [%s] %-58s %.3e\n", ok ? "PASS" : "FAIL", name.c_str(), value);
    std::fflush(stdout);
    return ok ? 0 : 1;
}

ProblemData make_pdata(int N, int Naxis) {
    ProblemData p;
    p.N = N;
    p.Naxis = Naxis;
    p.dt = 20e-6;
    const int Ntot = N * Naxis;
    p.X0.setZero(Ntot);
    p.inv_vec.setOnes(Ntot);
    p.set_vals.setOnes(Ntot);
    p.set_vals.array() *= NAN;
    for (int j = 0; j < Naxis; j++) {
        p.set_vals(j * N) = 0.0;
        p.set_vals(j * N + N - 1) = 0.0;
    }
    p.fixer.setOnes(Ntot);
    for (int j = 0; j < Naxis; j++) {
        p.fixer(j * N) = 0.0;
        p.fixer(j * N + N - 1) = 0.0;
    }
    return p;
}

} // namespace

int run_safe_signed_tests() {
    std::printf("\nOp_SAFE signed terms 1 and 3\n");
    int failures = 0;

    const int N = 64, Naxis = 3;
    ProblemData p = make_pdata(N, Naxis);

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

    return failures;
}
