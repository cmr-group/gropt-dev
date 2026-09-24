// Op_Slew::check and set_vals.
//
// check() skips a slew sample only when every waveform sample around it is pinned by set_vals, since a
// pinned-to-pinned step is not the solver's to fix. It used to look up set_vals with the slew vector's index
// on axis 0 only: with the combined norm (rot_variant false) a pinned axis switched the check off for the
// others (gx = 0 during a bipolar hid the bipolar's gy slew), and per axis every axis after the first read
// set_vals one sample early.

#include "Eigen/Dense"
#include "op_slew.hpp"
#include "problem_data.hpp"

#include <cmath>
#include <cstdio>
#include <string>

using namespace Gropt;

namespace {

int report(bool ok, const std::string &name, double value) {
    std::printf("  [%s] %-62s %.6g\n", ok ? "PASS" : "FAIL", name.c_str(), value);
    return ok ? 0 : 1;
}

ProblemData make_pdata(int N, int Naxis) {
    ProblemData p;
    p.N = N;
    p.Naxis = Naxis;
    p.dt = 10e-6;
    p.X0.setZero(N * Naxis);
    p.inv_vec.setOnes(N * Naxis);
    p.set_vals = Eigen::VectorXd::Constant(N * Naxis, NAN);
    p.fixer.setOnes(N * Naxis);
    return p;
}

// Waveform (ramps well under smax) with one step of `step` x smax on axis `ax` between samples i and i+1.
Eigen::VectorXd with_step(int N, int Naxis, double smax, double dt, int ax, int i, double step = 1.02) {
    Eigen::VectorXd X = Eigen::VectorXd::Zero(N * Naxis);
    for (int a = 0; a < Naxis; a++) {
        for (int j = 1; j < N; j++) X(a * N + j) = X(a * N + j - 1) + 0.3 * smax * dt * std::sin(0.2 * j + a);
    }
    const double d = step * smax * dt - (X(ax * N + i + 1) - X(ax * N + i));
    for (int j = i + 1; j < N; j++) X(ax * N + j) += d;
    return X;
}

// check() verdict for waveform X: 1 feasible, 0 not.
int verdict(ProblemData &p, bool rot_variant, double smax, const Eigen::VectorXd &X) {
    Op_Slew op(p, smax, rot_variant, 1.0);
    op.init();
    Eigen::VectorXd Xc = X, Ax(op.Ax_size);
    op.forward_op(Xc, Ax);
    op.check(Ax);
    return op.hist_feas.back();
}

} // namespace

int run_slew_check_tests() {
    std::printf("\n=== Op_Slew::check with set_vals ===\n");
    int fail = 0;
    const int N = 60;
    const double smax = 150.0, dt = 10e-6;

    // Combined norm: gx = 0 throughout and pinned on samples 0..19; a gy step at 10 -> 11 must be caught, and
    // the same waveform with a legal step passes.
    {
        ProblemData p = make_pdata(N, 2);
        for (int i = 0; i < 20; i++) p.set_vals(i) = 0.0;
        Eigen::VectorXd X = with_step(N, 2, smax, dt, 1, 10);
        Eigen::VectorXd Xok = with_step(N, 2, smax, dt, 1, 10, 0.9);
        X.head(N).setZero();
        Xok.head(N).setZero();
        fail += report(verdict(p, false, smax, Xok) == 1, "combined: legal waveform with gx pinned passes", 0);
        fail += report(verdict(p, false, smax, X) == 0, "combined: gx pinned does not hide a gy violation", 0);
    }
    // Per axis: axis 1 pinned on 30..39. A step 39 -> 40 leaves the pinned stretch (sample 40 free) and must be
    // caught; the old lookup read samples 37..39 and skipped it.
    {
        ProblemData p = make_pdata(N, 2);
        Eigen::VectorXd X = with_step(N, 2, smax, dt, 1, 39);
        for (int i = 30; i < 40; i++) p.set_vals(N + i) = X(N + i);
        fail += report(verdict(p, true, smax, X) == 0, "per axis: a step out of a pinned stretch is caught", 0);
    }
    // A step strictly inside a pinned stretch (samples i-1, i, i+1 all pinned) is skipped, on either axis.
    {
        ProblemData p = make_pdata(N, 2);
        Eigen::VectorXd X = with_step(N, 2, smax, dt, 1, 35);
        for (int i = 30; i < 40; i++) p.set_vals(N + i) = X(N + i);
        fail += report(verdict(p, true, smax, X) == 1, "per axis: a pinned-to-pinned step is skipped", 0);
        for (int i = 30; i < 40; i++) p.set_vals(i) = X(i);
        fail += report(verdict(p, false, smax, X) == 1, "combined: skipped only when every axis is pinned", 0);
    }
    // One axis: unchanged behaviour (skip inside a pinned stretch, catch next to a free sample).
    {
        ProblemData p = make_pdata(N, 1);
        Eigen::VectorXd X = with_step(N, 1, smax, dt, 0, 25);
        for (int i = 20; i < 30; i++) p.set_vals(i) = X(i);
        fail += report(verdict(p, true, smax, X) == 1, "one axis: pinned-to-pinned step skipped", 0);
        Eigen::VectorXd X2 = with_step(N, 1, smax, dt, 0, 29);
        for (int i = 20; i < 30; i++) p.set_vals(i) = X2(i);
        fail += report(verdict(p, true, smax, X2) == 0, "one axis: step next to a free sample caught", 0);
    }
    return fail;
}
