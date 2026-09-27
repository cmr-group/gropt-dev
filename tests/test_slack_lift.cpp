// The SAFE lifting: Op_SlackAbs (u >= |slew|) and Op_SAFE_Slack (term 2 reads u). Adjoints must be exact,
// and the lifted forward at u = |slew| must reproduce Op_SAFE (a reformulation, not a new constraint).

#include "op_diffbasin.hpp"
#include "op_gradient.hpp"
#include "op_moment.hpp"
#include "op_safe.hpp"
#include "op_safe_slack.hpp"
#include "op_slack_abs.hpp"
#include "test_util.hpp"

#include <random>

using namespace Gropt;
using namespace gropt_test;

namespace {

// Pinned ends, a 180 mid-way, and the slack block declared (as prepare() would) before any init().
ProblemData lifted_pdata(int N, int Naxis) {
    ProblemData p = make_pdata(N, Naxis, 20e-6, /*pin_ends=*/true);
    for (int j = 0; j < Naxis; j++) p.inv_vec.segment(j * N + N / 2, N - N / 2).setConstant(-1.0);
    p.add_aux(SAFE_SLACK_BLOCK, N * Naxis);
    p.fixer.conservativeResize(p.n_total());
    p.fixer.tail(N * Naxis).setOnes(); // auxiliary entries are free
    return p;
}

void random_fill(Eigen::VectorXd &v, std::mt19937 &rng) {
    std::normal_distribution<double> normal(0.0, 1.0);
    for (int i = 0; i < v.size(); i++) v(i) = normal(rng);
}

// Random waveform (times amp) with the slack block at |dx|, the tight value the coupling rows drive it to.
Eigen::VectorXd tight_point(const ProblemData &p, unsigned seed, double amp = 1.0) {
    std::mt19937 rng(seed);
    Eigen::VectorXd w(p.n_wave());
    random_fill(w, rng);
    Eigen::VectorXd X = Eigen::VectorXd::Zero(p.n_total());
    X.head(p.n_wave()) = amp * w;
    const int off = p.aux_offset(SAFE_SLACK_BLOCK);
    for (int j = 0; j < p.Naxis; j++) {
        for (int i = 0; i < p.N; i++) {
            const int k = j * p.N + i;
            X(off + k) = std::abs((i == 0) ? X(k) : (X(k) - X(k - 1)));
        }
    }
    return X;
}

// Raw forward/transpose see the operator's view of the primal: the waveform plus any aux block it declared.
int adjoint_raw(Operator &op, unsigned seed, const std::string &name) {
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
    return report(rel < 1e-12, name, rel);
}

} // namespace

int run_slack_lift_tests() {
    std::printf("\nSAFE lifting\n");
    int failures = 0;

    for (int Naxis : {1, 3}) {
        for (int N : {32, 257}) {
            ProblemData p = lifted_pdata(N, Naxis);
            const std::string tag = "N=" + std::to_string(N) + " Naxis=" + std::to_string(Naxis);

            // --- Op_SlackAbs: linear, so the raw adjoint must be exact ---
            Op_SlackAbs sa(p, SAFE_SLACK_BLOCK, 1.0);
            sa.init();
            failures += adjoint_raw(sa, 11u, "Op_SlackAbs adjoint, " + tag);

            // u >= |slew| holds with equality at the tight point, i.e. the smallest row is exactly 0
            {
                Eigen::VectorXd x = tight_point(p, 5u);
                Eigen::VectorXd ax = Eigen::VectorXd::Zero(sa.Ax_size);
                sa.forward(x, ax);
                failures += report(std::abs(ax.minCoeff()) < 1e-12, "tight slack sits ON its bound, " + tag,
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

                Eigen::VectorXd x = tight_point(p, 7u, 0.01); // waveform-sized amplitudes
                Eigen::VectorXd a_lift = Eigen::VectorXd::Zero(lifted.Ax_size);
                Eigen::VectorXd a_plain = Eigen::VectorXd::Zero(plain.Ax_size);
                lifted.forward(x, a_lift);
                plain.forward(x, a_plain);
                const double d = (a_lift - a_plain).cwiseAbs().maxCoeff();
                const double scale = a_plain.cwiseAbs().maxCoeff() + 1e-300;
                failures += report(d / scale < 1e-12,
                                   "lifted forward at u=|slew| equals Op_SAFE, " + tag, d / scale);

                // The CG sees the sign-frozen operator; its adjoint must be exact there.
                lifted.freeze_linearization(x);
                failures += adjoint_raw(lifted, 13u, "Op_SAFE_Slack frozen adjoint, " + tag);
                lifted.unfreeze_linearization();
            }
        }
    }

    // Ordinary operators get only the waveform slice, so a declared aux block must not change them.
    {
        ProblemData p = lifted_pdata(64, 1);
        Op_Gradient grad(p, 0.04, true, 1.0);
        grad.init();
        failures += adjoint_raw(grad, 21u, "Op_Gradient adjoint with an aux block");

        Op_Moment mom(p, 0.0, 0.0, 1e-4, "mT*ms/m", 0, 0, -1, 0, 1.0);
        mom.init();
        failures += adjoint_raw(mom, 22u, "Op_Moment adjoint with an aux block");

        Op_DiffBasin basin(p, 1e-3, 0.07, 0.04, 1.0, false);
        basin.init();
        failures += adjoint_raw(basin, 23u, "Op_DiffBasin adjoint with an aux block");
    }

    // The lifted op inherits the signed-terms prox (see test_safe_signed.cpp), so nothing is linearized.
    {
        ProblemData p = lifted_pdata(64, 3);
        Op_SAFE_Slack sgn(p, 1.0, 1.0, SAFE_SLACK_BLOCK);
        sgn.safe_params.set_demo_params();
        sgn.signed_terms13 = true;
        sgn.init();
        Eigen::VectorXd x = tight_point(p, 31u, 0.01);
        sgn.freeze_linearization(x);
        const double le = sgn.linearization_error(x);
        sgn.unfreeze_linearization();
        failures += report(le == 0.0, "lifted + signed terms: nothing is linearized at all", le);
    }

    return failures;
}
