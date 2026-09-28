// Shared helpers for the C++ tests. Header-only: every function is inline or a template.
#pragma once

#include "Eigen/Dense"
#include "problem_data.hpp"

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

namespace gropt_test {

constexpr double PI = 3.14159265358979323846; // MSVC does not define M_PI by default

// Print PASS/FAIL, the name and (unless NaN) a value; return 1 on failure.
inline int report(bool ok, const std::string &name, double value = NAN) {
    if (std::isnan(value)) std::printf("  [%s] %s\n", ok ? "PASS" : "FAIL", name.c_str());
    else std::printf("  [%s] %-62s %.6g\n", ok ? "PASS" : "FAIL", name.c_str(), value);
    std::fflush(stdout); // a crash in a later check must not hide the earlier ones
    return ok ? 0 : 1;
}

// Pass when within tol absolutely or relatively; tol = 0 means exactly equal.
inline int check_close(double actual, double expected, double tol, const std::string &name) {
    const double diff = std::abs(actual - expected), rel = diff / (std::abs(expected) + 1e-300);
    const bool ok = diff <= tol || rel <= tol;
    std::printf("  [%s] %-62s diff=%.2e  rel=%.2e\n", ok ? "PASS" : "FAIL", name.c_str(), diff, rel);
    std::fflush(stdout);
    return ok ? 0 : 1;
}

template <class F> bool throws_invalid(F &&f) {
    try {
        f();
    } catch (const std::invalid_argument &) {
        return true;
    }
    return false;
}

// X0 = 0, inv_vec = 1, all samples free; pin_ends fixes each axis's first and last sample to 0.
inline Gropt::ProblemData make_pdata(int N, int Naxis, double dt = 10e-6, bool pin_ends = false) {
    Gropt::ProblemData p;
    p.N = N;
    p.Naxis = Naxis;
    p.dt = dt;
    const int n = N * Naxis;
    p.X0.setZero(n);
    p.inv_vec.setOnes(n);
    p.set_vals = Eigen::VectorXd::Constant(n, NAN);
    p.fixer.setOnes(n);
    for (int j = 0; pin_ends && j < Naxis; j++) {
        p.set_vals(j * N) = p.set_vals(j * N + N - 1) = 0.0;
        p.fixer(j * N) = p.fixer(j * N + N - 1) = 0.0;
    }
    return p;
}

} // namespace gropt_test
