// Dot-product transpose tests for Op_* operators.
//
// For a linear operator A, <Ax, y> == <x, A^T y> for any x, y.
// Random x and y, check the residual is at machine precision.
//
// Three test sections:
//   A. Raw forward()/transpose() — operator-specific math, pure linear A
//   B. forward_op()/transpose_op() with apply_fixer=false — wrapper bookkeeping
//      (spec_norm, equilibration). Should also be a true transpose pair.
//   C. forward_op()/transpose_op() with default behavior (apply_fixer=true)
//      and a non-trivial fixer mask. The wrapper applies fixer to the
//      *output* of transpose but NOT to the input of forward, so the pair
//      is only a true transpose on the fixer-masked input subspace.
//      Two sub-cases per op:
//         C1. x masked by fixer first  -> expected pass
//         C2. x random/unmasked        -> expected fail (the asymmetry)

#include "Eigen/Dense"
#include "op_bvalue.hpp"
#include "op_gradient.hpp"
#include "op_slew.hpp"
#include "problem_data.hpp"

#include <cmath>
#include <cstdio>
#include <random>
#include <string>

using namespace Gropt;

namespace {

// Return 0 on pass, 1 on fail. Print PASS/FAIL either way.
int check_close(double actual, double expected, double tol, const std::string &name) {
    double diff = std::abs(actual - expected);
    double rel = diff / (std::abs(expected) + 1e-300);
    bool ok = diff < tol || rel < tol;
    std::printf("  [%s] %-60s diff=%.2e  rel=%.2e\n",
                ok ? "PASS" : "FAIL", name.c_str(), diff, rel);
    return ok ? 0 : 1;
}

// ProblemData with endpoints fixed (mimics vec_init_simple).
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

void random_fill(Eigen::VectorXd &v, std::mt19937 &rng) {
    std::normal_distribution<double> normal(0.0, 1.0);
    for (int i = 0; i < v.size(); i++) v(i) = normal(rng);
}

// --- A: raw forward/transpose dot-product test ---
int dot_product_raw(Operator &op, int Ntot, int Ax_size, unsigned seed,
                    const std::string &name, double tol = 1e-12) {
    std::mt19937 rng(seed);
    Eigen::VectorXd x(Ntot), y(Ax_size);
    random_fill(x, rng);
    random_fill(y, rng);

    Eigen::VectorXd Ax(Ax_size), Aty(Ntot);
    op.forward(x, Ax);
    op.transpose(y, Aty);
    return check_close(Ax.dot(y), x.dot(Aty), tol, name);
}

// --- B: forward_op/transpose_op with apply_fixer=false ---
int dot_product_op_no_fixer(Operator &op, int Ntot, int Ax_size, unsigned seed,
                            const std::string &name, double tol = 1e-12) {
    std::mt19937 rng(seed);
    Eigen::VectorXd x(Ntot), y(Ax_size);
    random_fill(x, rng);
    random_fill(y, rng);

    Eigen::VectorXd Ax(Ax_size), Aty(Ntot);
    op.forward_op(x, Ax);
    op.transpose_op(y, Aty, /*apply_fixer=*/false);
    return check_close(Ax.dot(y), x.dot(Aty), tol, name);
}

// --- C: forward_op/transpose_op with default (apply_fixer=true) ---
// `mask_input` controls whether x is multiplied by fixer before the test.
int dot_product_op_with_fixer(Operator &op, int Ntot, int Ax_size, unsigned seed,
                              bool mask_input, const std::string &name,
                              double tol = 1e-12) {
    std::mt19937 rng(seed);
    Eigen::VectorXd x(Ntot), y(Ax_size);
    random_fill(x, rng);
    random_fill(y, rng);

    if (mask_input) {
        x.array() *= op.pdata->fixer.array();
    }

    Eigen::VectorXd Ax(Ax_size), Aty(Ntot);
    op.forward_op(x, Ax);
    op.transpose_op(y, Aty); // default apply_fixer=true
    return check_close(Ax.dot(y), x.dot(Aty), tol, name);
}

// Helper: build, init, and dispatch raw dot-product tests over a sweep.
template <typename OpFactory>
int sweep_raw(const std::string &op_name, OpFactory make_op) {
    int failures = 0;
    const int sizes[] = {16, 64, 257};
    const int axes[] = {1, 3};
    for (int N : sizes) {
        for (int Naxis : axes) {
            ProblemData pdata = make_pdata(N, Naxis);
            auto op = make_op(pdata);
            op->init();
            for (unsigned seed = 1; seed <= 3; seed++) {
                char label[96];
                std::snprintf(label, sizeof(label), "%s raw  N=%3d Naxis=%d seed=%u",
                              op_name.c_str(), N, Naxis, seed);
                failures += dot_product_raw(*op, op->Ntot, op->Ax_size, seed, label);
            }
        }
    }
    return failures;
}

// Helper: smaller sweep for the _op variants (just one shape, two seeds).
template <typename OpFactory>
int sweep_op_variants(const std::string &op_name, OpFactory make_op) {
    int failures = 0;
    const int N = 64;
    const int Naxis = 3;
    ProblemData pdata = make_pdata(N, Naxis);
    auto op = make_op(pdata);
    op->init();

    for (unsigned seed = 1; seed <= 2; seed++) {
        // B: wrapper without fixer — should pass
        {
            char label[128];
            std::snprintf(label, sizeof(label), "%s _op  no-fixer        N=%d Naxis=%d seed=%u",
                          op_name.c_str(), N, Naxis, seed);
            failures += dot_product_op_no_fixer(*op, op->Ntot, op->Ax_size, seed, label);
        }
        // C1: wrapper with fixer, x pre-masked — should pass
        {
            char label[128];
            std::snprintf(label, sizeof(label), "%s _op  fixer + masked x N=%d Naxis=%d seed=%u",
                          op_name.c_str(), N, Naxis, seed);
            failures += dot_product_op_with_fixer(*op, op->Ntot, op->Ax_size, seed,
                                                  /*mask_input=*/true, label);
        }
        // C2: wrapper with fixer, x random — INTENTIONALLY FAILS, kept as
        // documentation of the wrapper's asymmetry (transpose_op applies fixer
        // to its output, forward_op does not). Re-enable to confirm the gap.
        // {
        //     char label[128];
        //     std::snprintf(label, sizeof(label), "%s _op  fixer + raw x    N=%d Naxis=%d seed=%u",
        //                   op_name.c_str(), N, Naxis, seed);
        //     failures += dot_product_op_with_fixer(*op, op->Ntot, op->Ax_size, seed,
        //                                           /*mask_input=*/false, label);
        // }
    }
    return failures;
}

} // namespace

int run_op_transpose_tests() {
    int failures = 0;

    auto make_bvalue = [](ProblemData &p) {
        return std::make_unique<Op_BValue>(p, /*bval_target=*/100.0, /*bval_tol0=*/10.0,
                                           /*start_idx0=*/-1, /*stop_idx0=*/-1,
                                           /*weight_mod=*/1.0,
                                           BVALUE_MODE_MINVAL,
                                           /*max_scale=*/1.01);
    };

    std::printf("\n=== A. Raw forward/transpose ===\n");
    failures += sweep_raw("Op_Gradient", [](ProblemData &p) {
        return std::make_unique<Op_Gradient>(p, /*gmax=*/0.04, /*rot_variant=*/true, /*weight_mod=*/1.0);
    });
    failures += sweep_raw("Op_Slew", [](ProblemData &p) {
        return std::make_unique<Op_Slew>(p, /*smax=*/200.0, /*rot_variant=*/true, /*weight_mod=*/1.0);
    });
    failures += sweep_raw("Op_BValue", make_bvalue);

    std::printf("\n=== B+C. forward_op / transpose_op variants ===\n");
    failures += sweep_op_variants("Op_Gradient", [](ProblemData &p) {
        return std::make_unique<Op_Gradient>(p, /*gmax=*/0.04, /*rot_variant=*/true, /*weight_mod=*/1.0);
    });
    failures += sweep_op_variants("Op_Slew", [](ProblemData &p) {
        return std::make_unique<Op_Slew>(p, /*smax=*/200.0, /*rot_variant=*/true, /*weight_mod=*/1.0);
    });
    failures += sweep_op_variants("Op_BValue", make_bvalue);

    return failures;
}
