// Gradient-direction tests for Operator::add_obj.
//
// add_obj implements the LHS contribution of a quadratic objective
// f(X) = (1/2) ||A X||^2 — its gradient is A^T A X, so add_obj should add
// obj_weight * A^T A X to the linear-system output vector.
//
// We compute the expected vector explicitly by chaining forward_op +
// transpose_op, then compare against what add_obj actually produces.
//
// Two operators:
//   - Op_Identity (A = I): A^T A X == X == A X. add_obj is correct
//     regardless of which intermediate vector it uses, so this is a sanity
//     check that the test scaffolding works.
//   - Op_BValue (A is the cumulative-integration operator): A^T A X != A X.
//     This is the discriminating case — if add_obj uses the wrong vector,
//     this is where it shows up.

#include "Eigen/Dense"
#include "op_bvalue.hpp"
#include "op_identity.hpp"
#include "problem_data.hpp"

#include <cmath>
#include <cstdio>
#include <memory>
#include <random>
#include <string>

using namespace Gropt;

namespace {

// Compare two vectors element-wise. Also reports cosine similarity, which is
// the most readable signal for "wrong direction" (the failure mode we expect
// if add_obj uses A x instead of A^T A x).
int check_close_vec(const Eigen::VectorXd &actual, const Eigen::VectorXd &expected,
                    double tol, const std::string &name) {
    Eigen::VectorXd diff = actual - expected;
    double linf = diff.cwiseAbs().maxCoeff();
    double rel = linf / (expected.cwiseAbs().maxCoeff() + 1e-300);
    bool ok = linf < tol || rel < tol;

    double na = actual.norm();
    double ne = expected.norm();
    double cos_theta = actual.dot(expected) / ((na * ne) + 1e-300);

    std::printf("  [%s] %-56s linf=%.2e  rel=%.2e  cos=%+.6f\n",
                ok ? "PASS" : "FAIL", name.c_str(), linf, rel, cos_theta);
    return ok ? 0 : 1;
}

// Minimal pdata: no endpoint fixing, identity inv_vec. Keeps the test focused
// on the add_obj math instead of fixer/equilibration interactions.
ProblemData make_pdata(int N, int Naxis) {
    ProblemData p;
    p.N = N;
    p.Naxis = Naxis;
    p.dt = 10e-6;

    int Ntot = N * Naxis;
    p.X0.setZero(Ntot);
    p.inv_vec.setOnes(Ntot);
    p.set_vals.setOnes(Ntot);
    p.set_vals.array() *= NAN;
    p.fixer.setOnes(Ntot);
    return p;
}

void random_fill(Eigen::VectorXd &v, std::mt19937 &rng) {
    std::normal_distribution<double> normal(0.0, 1.0);
    for (int i = 0; i < v.size(); i++) v(i) = normal(rng);
}

// Expected LHS contribution: obj_weight * (forward_op then transpose_op)(X).
// This is the gradient of (1/2) ||A X||^2 scaled by obj_weight, expressed
// using exactly the same wrapper calls add_obj uses internally.
Eigen::VectorXd expected_grad(Operator &op, const Eigen::VectorXd &X) {
    Eigen::VectorXd X_in = X; // forward_op takes a non-const ref
    Eigen::VectorXd Ax(op.Ax_size);
    Ax.setZero();
    op.forward_op(X_in, Ax);

    Eigen::VectorXd AtAx(op.Ntot);
    AtAx.setZero();
    op.transpose_op(Ax, AtAx);

    return op.obj_weight * AtAx;
}

// Actual LHS contribution: zero out, call add_obj, return what got added.
Eigen::VectorXd actual_add_obj(Operator &op, const Eigen::VectorXd &X) {
    Eigen::VectorXd X_in = X;
    Eigen::VectorXd out(op.Ntot);
    out.setZero();
    op.add_obj(X_in, out);
    return out;
}

template <typename OpFactory>
int sweep(const std::string &op_name, OpFactory make_op, double tol = 1e-12) {
    int failures = 0;
    const int sizes[] = {16, 64, 257};
    const int axes[] = {1, 3};

    for (int N : sizes) {
        for (int Naxis : axes) {
            ProblemData pdata = make_pdata(N, Naxis);
            auto op = make_op(pdata);
            op->init();

            std::mt19937 rng(42u + static_cast<unsigned>(N) + 7u * static_cast<unsigned>(Naxis));
            for (unsigned seed = 1; seed <= 2; seed++) {
                Eigen::VectorXd X(op->Ntot);
                random_fill(X, rng);

                Eigen::VectorXd expected = expected_grad(*op, X);
                Eigen::VectorXd actual = actual_add_obj(*op, X);

                char label[96];
                std::snprintf(label, sizeof(label), "%s add_obj  N=%3d Naxis=%d seed=%u",
                              op_name.c_str(), N, Naxis, seed);
                failures += check_close_vec(actual, expected, tol, label);
            }
        }
    }
    return failures;
}

} // namespace

int run_add_obj_tests() {
    int failures = 0;

    std::printf("\n=== add_obj should produce obj_weight * A^T A X ===\n");
    std::printf("(Op_Identity has A = I, so any add_obj that touches forward+transpose\n");
    std::printf(" will pass. Op_BValue has A != I; if add_obj returns A x rather than\n");
    std::printf(" A^T A x, the BValue cases fail and cos(actual,expected) drops below 1.)\n");

    failures += sweep("Op_Identity", [](ProblemData &p) {
        return std::make_unique<Op_Identity>(p, /*weight_mod=*/1.0);
    });
    failures += sweep("Op_BValue", [](ProblemData &p) {
        return std::make_unique<Op_BValue>(p, /*bval_target=*/100.0,
                                           /*bval_tol0=*/10.0,
                                           /*start_idx0=*/-1, /*stop_idx0=*/-1,
                                           /*weight_mod=*/1.0,
                                           BVALUE_MODE_MINVAL,
                                           /*max_scale=*/1.01);
    });

    return failures;
}
