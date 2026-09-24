// The SAFE lifting: Op_SlackAbs (u >= |slew|) and Op_SAFE_Slack (term 2 reads u).
//
// Two things have to hold. The adjoints must be exact, or the CG solves a different problem than the one
// the dual update assumes. And the lifted forward evaluated at u = |slew| must reproduce Op_SAFE exactly:
// that is what makes the lifting a reformulation rather than a different constraint.

#include "Eigen/Dense"
#include "op_diffbasin.hpp"
#include "op_gradient.hpp"
#include "op_moment.hpp"
#include "op_safe.hpp"
#include "op_safe_slack.hpp"
#include "op_slack_abs.hpp"
#include "problem_data.hpp"

#include <cmath>
#include <cstdio>
#include <memory>
#include <random>
#include <string>

using namespace Gropt;

namespace {

int check(bool ok, const std::string &name, double value) {
    std::printf("  [%s] %-60s %.3e\n", ok ? "PASS" : "FAIL", name.c_str(), value);
    std::fflush(stdout); // a crash in a later check must not hide the earlier ones
    return ok ? 0 : 1;
}

ProblemData make_pdata(int N, int Naxis) {
    ProblemData p;
    p.N = N;
    p.Naxis = Naxis;
    p.dt = 20e-6;
    p.add_aux(SAFE_SLACK_BLOCK, N * Naxis); // as GroptParams::prepare() would, before any init()

    p.X0.setZero(p.n_wave());
    p.inv_vec.setOnes(p.n_wave());
    for (int j = 0; j < Naxis; j++) { // a 180 in the middle, as diffusion has
        for (int i = N / 2; i < N; i++) p.inv_vec(j * N + i) = -1.0;
    }
    p.set_vals.setOnes(p.n_wave());
    p.set_vals.array() *= NAN;
    for (int j = 0; j < Naxis; j++) {
        p.set_vals(j * N) = 0.0;
        p.set_vals(j * N + N - 1) = 0.0;
    }
    p.fixer.setOnes(p.n_total()); // auxiliary entries are free
    for (int j = 0; j < Naxis; j++) {
        p.fixer(j * N) = 0.0;
        p.fixer(j * N + N - 1) = 0.0;
    }
    return p;
}

void random_fill(Eigen::VectorXd &v, std::mt19937 &rng) {
    std::normal_distribution<double> normal(0.0, 1.0);
    for (int i = 0; i < v.size(); i++) v(i) = normal(rng);
}

// Fill the slack block with |dx|, the tight value the coupling rows drive it to (waveform units).
void set_tight_slack(Eigen::VectorXd &X, const ProblemData &p) {
    const int off = p.aux_offset(SAFE_SLACK_BLOCK);
    for (int j = 0; j < p.Naxis; j++) {
        for (int i = 0; i < p.N; i++) {
            const int k = j * p.N + i;
            const double d = (i == 0) ? X(k) : (X(k) - X(k - 1));
            X(off + k) = std::abs(d);
        }
    }
}

// The raw forward/transpose of an operator work over ITS view of the primal: the waveform, unless it
// declared an auxiliary block.
int adjoint_raw(Operator &op, const ProblemData &p, unsigned seed, const std::string &name) {
    std::mt19937 rng(seed);
    Eigen::VectorXd x(op.n_primal()), y(op.Ax_size);
    random_fill(x, rng);
    random_fill(y, rng);
    Eigen::VectorXd Ax = Eigen::VectorXd::Zero(op.Ax_size);
    Eigen::VectorXd Aty = Eigen::VectorXd::Zero(op.n_primal());
    op.forward(x, Ax);
    op.transpose(y, Aty);
    const double lhs = Ax.dot(y), rhs = x.dot(Aty);
    const double rel = std::abs(lhs - rhs) / (std::abs(lhs) + 1e-300);
    return check(rel < 1e-12, name, rel);
}

} // namespace

int run_slack_lift_tests() {
    std::printf("\nSAFE lifting\n");
    int failures = 0;

    for (int Naxis : {1, 3}) {
        for (int N : {32, 257}) {
            ProblemData p = make_pdata(N, Naxis);
            const std::string tag = "N=" + std::to_string(N) + " Naxis=" + std::to_string(Naxis);

            // --- Op_SlackAbs: linear, so the raw adjoint must be exact ---
            Op_SlackAbs sa(p, SAFE_SLACK_BLOCK, 1.0);
            sa.init();
            failures += adjoint_raw(sa, p, 11u, "Op_SlackAbs adjoint, " + tag);

            // u >= |slew| must hold with equality at the tight point, i.e. one row is exactly 0
            {
                std::mt19937 rng(5u);
                Eigen::VectorXd x = Eigen::VectorXd::Zero(p.n_total());
                Eigen::VectorXd w(p.n_wave());
                random_fill(w, rng);
                x.head(p.n_wave()) = w;
                set_tight_slack(x, p);
                Eigen::VectorXd ax = Eigen::VectorXd::Zero(sa.Ax_size);
                sa.forward(x, ax);
                failures += check(ax.minCoeff() >= -1e-12, "tight slack is feasible, " + tag,
                                  ax.minCoeff());
                failures += check(ax.minCoeff() < 1e-12, "tight slack is ON a bound, " + tag,
                                  ax.minCoeff());
            }

            // --- Op_SAFE_Slack at u = |slew| must equal Op_SAFE exactly ---
            {
                Op_SAFE_Slack lifted(p, 1.0, 1.0, SAFE_SLACK_BLOCK);
                lifted.safe_params.set_demo_params();
                lifted.init();
                Op_SAFE plain(p, 1.0, 1.0);
                plain.safe_params.set_demo_params();
                plain.init();

                std::mt19937 rng(7u);
                Eigen::VectorXd x = Eigen::VectorXd::Zero(p.n_total());
                Eigen::VectorXd w(p.n_wave());
                random_fill(w, rng);
                w *= 0.01; // waveform-sized amplitudes
                x.head(p.n_wave()) = w;
                set_tight_slack(x, p);

                Eigen::VectorXd a_lift = Eigen::VectorXd::Zero(lifted.Ax_size);
                Eigen::VectorXd a_plain = Eigen::VectorXd::Zero(plain.Ax_size);
                lifted.forward(x, a_lift);
                plain.forward(x, a_plain);
                const double d = (a_lift - a_plain).cwiseAbs().maxCoeff();
                const double scale = a_plain.cwiseAbs().maxCoeff() + 1e-300;
                failures += check(d / scale < 1e-12,
                                  "lifted forward at u=|slew| equals Op_SAFE, " + tag, d / scale);

                // The CG sees the sign-frozen operator; its adjoint must be exact there.
                lifted.freeze_linearization(x);
                failures += adjoint_raw(lifted, p, 13u, "Op_SAFE_Slack frozen adjoint, " + tag);
                lifted.unfreeze_linearization();
            }
        }
    }


    // Ordinary operators must be COMPLETELY unaffected by a declared auxiliary block: they are handed
    // the waveform slice, so their whole-vector expressions (out = X, A * X, a.dot(X)) still see exactly
    // the vector they always did. These three are the ones that would break if they ever saw the longer
    // primal, so they are the sharpest check that the plumbing keeps it away from them.
    {
        ProblemData p = make_pdata(64, 1);
        Op_Gradient grad(p, 0.04, true, 1.0);
        grad.init();
        failures += adjoint_raw(grad, p, 21u, "Op_Gradient adjoint with an aux block");

        Op_Moment mom(p, 0.0, 0.0, 1e-4, "mT*ms/m", 0, 0, -1, 0, 1.0);
        mom.init();
        failures += adjoint_raw(mom, p, 22u, "Op_Moment adjoint with an aux block");

        Op_DiffBasin basin(p, 1e-3, 0.07, 0.04, 1.0, false);
        basin.init();
        failures += adjoint_raw(basin, p, 23u, "Op_DiffBasin adjoint with an aux block");
    }

    // The signed-terms projection itself is covered by tests/test_safe_signed.cpp; what matters
    // here is only that the lifted operator inherits it, which the lin_err check below pins.
    {
        ProblemData p = make_pdata(64, 3);
        Op_SAFE_Slack sgn(p, 1.0, 1.0, SAFE_SLACK_BLOCK);
        sgn.safe_params.set_demo_params();
        sgn.signed_terms13 = true;
        sgn.init();
        std::mt19937 rng(31u);
        Eigen::VectorXd x = Eigen::VectorXd::Zero(p.n_total());
        Eigen::VectorXd w(p.n_wave());
        random_fill(w, rng);
        w *= 0.01;
        x.head(p.n_wave()) = w;
        set_tight_slack(x, p);
        sgn.freeze_linearization(x);
        const double le = sgn.linearization_error(x);
        sgn.unfreeze_linearization();
        failures += check(le == 0.0,
                          "lifted + signed terms: nothing is linearized at all", le);
    }

    return failures;
}
