// demo_diffusion_simple.cpp -- diffusion waveforms from C++, start to finish: PNS + CNS + concomitant-field
// balance, M0..M2 nulled, b-value maximized. Solved on a coarse raster with a tuned recipe, then moved to the
// 10 us hardware raster and re-checked there.
//
// main() shows the two uses: solve_at(TE) for one waveform at a fixed TE, and te_search(b) for the shortest TE
// that reaches a b-value. Each writes its 10 us waveform to a text file.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <fstream>
#include <stdexcept>

#include "spdlog/spdlog.h"

#include "diffusion_recipe.hpp" // Recipe, load_recipe, set_recipe_x0, configure_solver
#include "gropt_params.hpp"     // GroptParams: the problem (timing + constraints)
#include "gropt_utils.hpp"      // get_SAFE_eigen
#include "op_bvalue.hpp"        // BVALUE_MODE_MINVAL
#include "op_safe_slack.hpp"    // SAFE_SLACK_BLOCK
#include "raster.hpp"           // raster_ratio, resample_waveform, resample_inv_vec
#include "solver_groptsdmm.hpp" // SolverGroptSDMM: the solver

using namespace Gropt;

namespace {

//  problem (everything but TE)
// -------------------------------------------------------------- 
const double T_90 = 3e-3, T_180 = 5e-3, T_readout = 16e-3;  // [s]; T_readout = readout time before TE
const double dt = 200e-6;      // solve raster [s]; the recipes are tuned at 100 to 400 us
const double dt_hw = 10e-6;    // hardware raster [s]; dt must be an integer multiple of it
const double gmax = 0.19;      // [T/m]
const double smax = 200.0;     // [T/m/s]
const double pns_lim = 0.8;    // fraction of the PNS threshold
const double cns_lim = 0.8;    // fraction of the cardiac threshold
const int MMT = 2;             // null M0..M_MMT
const bool do_concomitant = true; // balance concomitant field at TE
const bool safe_all_axes = false; // SAFE on x, y and z (safe on any axis); false = x only

// This section is just copyiing the random parameters that I use in python, in reality it needs to grab the scanner SAFE values.
// --------------------------------------------------------------
struct SafeTable {
    Eigen::VectorXd tau1{3}, tau2{3}, tau3{3}, a1{3}, a2{3}, a3{3}, stim_limit{3}, g_scale{3};
};

// Peripheral nerve stimulation model
SafeTable pns_table() {
    SafeTable s;
    s.tau1 << 0.86e-3, 0.93e-3, 0.78e-3;
    s.tau2 << 11.0e-3, 13.0e-3, 10.0e-3;
    s.tau3 << 0.27e-3, 0.23e-3, 0.25e-3;
    s.a1 << 0.28, 0.24, 0.29;
    s.a2 << 0.54, 0.42, 0.60;
    s.a3 << 0.18, 0.34, 0.11;
    s.stim_limit << 29.0, 27.4, 38.5;
    s.g_scale << 0.34, 0.34, 0.29;
    return s;
}

// Cardiac stimulation model
SafeTable cns_table() {
    SafeTable s;
    s.tau1 << 2.7e-3, 3.0e-3, 2.3e-3;
    s.tau2 << 1.4e-3, 1.5e-3, 1.2e-3;
    s.tau3 << 1.1e-3, 1.5e-3, 1.2e-3;
    s.a1 << 0.67, 0.79, 0.78;
    s.a2 << 0.33, 0.21, 0.22;
    s.a3 << 0.00, 0.00, 0.00;
    s.stim_limit << 14.3, 14.9, 18.1;
    s.g_scale << 0.34, 0.30, 0.32;
    return s;
}

const SafeTable PNS = pns_table(), CNS = cns_table();

// SAFE table axes the single-axis waveform is tested on: 0 (x), or 0, 1, 2 (x, y, z)
const int n_safe_axes = safe_all_axes ? 3 : 1;

// One SAFE limit, added once per tested axis. Passing on x, y and z also makes it safe along any unit diffusion
// direction, since the scanner combines axes as a root-sum-square.
void add_SAFE(GroptParams &gp, double lim, const SafeTable &s, double weight) {
    for (int ax = 0; ax < n_safe_axes; ++ax)
        gp.add_SAFE(lim, s.tau1, s.tau2, s.tau3, s.a1, s.a2, s.a3, s.stim_limit, s.g_scale, ax, weight);
}

// Peak SAFE response of g over the tested axes, as a fraction of the model's threshold, with the exact
// (raster-invariant) filter. This is a check at the end, not part of the optimization.
double peak_SAFE(const Eigen::VectorXd &g, double dt, const SafeTable &s) {
    double peak = 0.0;
    for (int ax = 0; ax < n_safe_axes; ++ax)
        peak = std::max(peak, get_SAFE_eigen(g, 1, dt, true, ax, s.tau1, s.tau2, s.tau3, s.a1, s.a2, s.a3,
                                             s.stim_limit, s.g_scale, true).cwiseAbs().maxCoeff());
    return peak;
}


// Result at one TE: the solve on the dt raster, the same waveform on the 10 us raster, and its checks there
// --------------------------------------------------------------
struct Result {
    double TE = 0.0;
    SolveResult solve;             // the solve on the dt raster
    Eigen::VectorXd g, inv;        // on the dt_hw raster: G [T/m], encoding sign
    double r_gmax = 0, r_slew = 0, r_pns = 0, r_cns = 0;   // peak / limit
    double conc = 0;               // |E_pre / E_post - 1|
    Eigen::VectorXd M;             // M0..M_MMT [mT/m ms^(k+1)]
    double b = 0;                  // [s/mm^2]
    bool ok = false;               // converged, every limit within 2 %, concomitant within its tol
};

// Build, solve and re-check one waveform at this TE. Each R.w_* is the starting weight the recipe tuned for
// that constraint; the limits themselves come from the problem above.
Result solve_at(double TE, const Recipe &R) {
    Result res;
    res.TE = TE;

    // ---- constraints ---------------------------------------------------------------------------------------
    GroptParams gp;
    gp.diff_init(dt, TE, T_90, T_180, T_readout);   // RF gaps, the 180 sign flip, zero at both ends
    gp.add_gmax(gmax, true, R.w_gmax);              // |G| <= gmax
    gp.add_smax(smax, true, R.w_smax);              // |dG/dt| <= smax

    // Set moment constraints
    // moment_project = true (recipe) enforces them exactly by projection.
    for (int m = 0; m <= MMT; ++m)
        gp.add_moment(m, 0.0, 1e-5, "mT*ms/m", 0, 0, -1, 0, R.w_moment, R.moment_project);

    // add_SAFE copies these into each SAFE operator, so set them first
    gp.safe_alpha_exact = true;    // the raster-invariant SAFE filter: off by default in C++, needed for 10 us
    
    // Everything else SAFE wise is tuned by the recipe, no real reason to change anything here:
    gp.safe_lifted = R.safe_lifted;       // the recipe's SAFE formulation (lifted = an extra |slew| variable)
    gp.safe_signed13 = R.safe_signed13;
    gp.safe_eps = R.safe_eps;
    if (R.safe_lifted) gp.ensure_slack_abs(SAFE_SLACK_BLOCK, R.w_slack);   // that variable, before add_SAFE
    add_SAFE(gp, pns_lim, PNS, R.w_pns);  // PNS <= pns_lim
    add_SAFE(gp, cns_lim, CNS, R.w_cns);  // CNS <= cns_lim

    // Whether to add concomittant constraint
    if (do_concomitant) {
        gp.add_concomitant(0, true, R.w_concomitant, 0.1, 1.0, true);
    }

    // Maximize b over the whole waveform (samples 0 to -1 = end). add_bvalue also builds the b-value
    // constraint; as an objective it ignores the constraint's target, tol, mode and max_scale (placeholders here).
    gp.add_bvalue(/*target*/ 100.0, /*tol*/ 1.0, 0, -1, R.bval_obj_weight, /*mode*/ BVALUE_MODE_MINVAL,
                  /*max_scale*/ 1.01, /*as_objective*/ true, /*linearize*/ true);
    gp.normalize_obj = true;       // b pull scaled to bval_obj_weight, as the recipes were tuned

    set_recipe_x0(gp, R, MMT);     // the recipe's starting waveform
    gp.prepare();                  // finalize the setup: call after every add_*

    // --------------------------------------------------------------------------------------------------------
    // ---- solve ---------------------------------------------------------------------------------------------
    SolverGroptSDMM solver;
    configure_solver(R, solver);   // the recipe's solver settings (iterations, reweighting, ...)
    res.solve = solver.solve(gp);  // .X = waveform on the dt raster [T/m]

    // ---- rasmple to hardware raster ------------------------------------------------------------------------
    // The scanner plays the line through the samples, so linear resampling onto dt_hw is the same waveform.
    const int ratio = raster_ratio(dt, dt_hw);                   // 20; throws if not an integer
    res.g = resample_waveform(res.solve.X, 1, ratio);            // 1 axis
    res.inv = resample_inv_vec(gp.pdata.inv_vec, 1, ratio);      // +1 before the 180, -1 after
    const Eigen::VectorXd &g = res.g, &inv = res.inv;
    const int n = static_cast<int>(g.size());

    //  All below is just checking, we are already done and g is out output.
    // ---- re-check every limit on the 10 us waveform, independently of the solver ---------------------------
    res.r_gmax = g.cwiseAbs().maxCoeff() / gmax;
    res.r_slew = (g.tail(n - 1) - g.head(n - 1)).cwiseAbs().maxCoeff() / dt_hw / smax;
    res.r_pns = peak_SAFE(g, dt_hw, PNS) / pns_lim;
    res.r_cns = peak_SAFE(g, dt_hw, CNS) / cns_lim;

    // Concomitant: energy before vs after the 180 (the RF gap is zero, so it falls in neither)
    double e_pre = 0.0, e_post = 0.0;
    for (int i = 0; i + 1 < n; ++i) {
        const double e = dt_hw / 3.0 * (g(i) * g(i) + g(i) * g(i + 1) + g(i + 1) * g(i + 1));  // exact per piece
        if (inv(i) > 0 && inv(i + 1) > 0) e_pre += e;
        if (inv(i) < 0 && inv(i + 1) < 0) e_post += e;
    }
    res.conc = std::abs(e_pre / e_post - 1.0);

    // Moments Mk = integral(inv * g * t^k); the rectangle rule is exact here (M0 = 0)
    res.M = Eigen::VectorXd::Zero(MMT + 1);
    for (int k = 0; k <= MMT; ++k) {
        for (int i = 0; i < n; ++i) res.M(k) += inv(i) * g(i) * std::pow(i * dt_hw, k) * dt_hw;
        res.M(k) *= 1e3 * std::pow(1e3, k + 1);   // T/m s^(k+1) -> mT/m ms^(k+1)
    }

    // b = gamma^2 * integral(F^2), F = integral(inv * g)
    double F = 0.0;
    for (int i = 0; i < n; ++i) {
        F += inv(i) * g(i) * dt_hw;
        res.b += F * F * dt_hw;
    }
    res.b *= 267.5221900e6 * 267.5221900e6 * 1e-6;   // gamma [rad/s/T]; s/m^2 -> s/mm^2

    // 2 % headroom on the limits, as Python's verify()
    res.ok = res.solve.converged && res.r_gmax <= 1.02 && res.r_slew <= 1.02 && res.r_pns <= 1.02 &&
             res.r_cns <= 1.02 && (!do_concomitant || res.conc <= 0.1);
    return res;
}

// Shortest TE that passes every check and reaches b_target on the 10 us raster. TE is quantized to the
// dt grid: every TE in [T_readout + m dt, T_readout + (m + 1) dt) gives the same N = m + 1 samples, so search
// over m and solve each bin at its middle. Longer TE gives more b, so bisect between an infeasible lo and a
// feasible hi (the solver is local, so b(TE) is only nearly monotone; the answer can be a bin or so off).
Result te_search(double b_target, const Recipe &R, double TE_max = 150e-3) {
    // Bin m -> the TE solved for it (the middle, so rounding never moves it to another bin)
    auto TE_of = [](int m) { return T_readout + (m + 0.5) * dt; };

    // Solve bin m; feasible = every check passes and b reaches the target
    auto evaluate = [&](int m, bool &feasible) {
        Result res = solve_at(TE_of(m), R);
        feasible = res.ok && res.b >= b_target;
        std::printf("  TE %6.2f ms: b %7.1f  %s\n", res.TE * 1e3, res.b, feasible ? "feasible" : "infeasible");
        return res;
    };

    // Shortest TE that fits the RF and readout. diff_init cannot build anything shorter, so start there.
    const double TE_min = std::max(T_90 + T_180 + T_readout, 2.0 * (T_180 + T_readout));
    int lo = static_cast<int>(std::ceil((TE_min - T_readout) / dt - 0.5));    // first bin with TE >= TE_min
    int hi = static_cast<int>(std::floor((TE_max - T_readout) / dt - 0.5));   // last bin with TE <= TE_max

    // The upper end must reach the target, else there is nothing to search
    bool feasible = false;
    Result best = evaluate(hi, feasible);
    if (!feasible) throw std::runtime_error("b target not reached at TE_max; raise TE_max");

    // A low target may already be reached at the shortest TE
    Result shortest = evaluate(lo, feasible);
    if (feasible) return shortest;

    // Halve the bracket until lo and hi are neighbouring bins; ~log2(number of bins) solves
    while (hi - lo > 1) {                 // lo infeasible, hi feasible
        const int mid = (lo + hi) / 2;
        Result res = evaluate(mid, feasible);
        if (feasible) { hi = mid; best = std::move(res); }   // keep the shortest feasible result so far
        else lo = mid;
    }
    return best;
}

// Print the solve and the 10 us re-check
void report(const Result &res) {
    std::printf("solve at %.0f us: TE = %.2f ms, N = %d, b = %.1f s/mm^2, converged = %s, %d iterations\n",
                dt * 1e6, res.TE * 1e3, static_cast<int>(res.solve.X.size()), res.solve.bvalue,
                res.solve.converged ? "yes" : "no", res.solve.n_iter);
    std::printf("\non the %.0f us raster (N = %d):\n", dt_hw * 1e6, static_cast<int>(res.g.size()));
    std::printf("  gmax   %.3f of the limit\n", res.r_gmax);
    std::printf("  slew   %.3f\n", res.r_slew);
    const char *axes = safe_all_axes ? "worst of x, y, z" : "x";
    std::printf("  PNS    %.3f  (%s)\n", res.r_pns, axes);
    std::printf("  CNS    %.3f  (%s)\n", res.r_cns, axes);
    std::printf("  concomitant  |E_pre / E_post - 1| = %.1e%s\n", res.conc,
                do_concomitant ? "" : "  (not constrained)");
    for (int k = 0; k <= MMT; ++k) std::printf("  M%d     %+.1e mT/m ms^%d\n", k, res.M(k), k + 1);
    std::printf("  b      %.1f s/mm^2\n", res.b);
}

// Save the 10 us waveform, one sample per line: G [T/m], encoding sign
void save(const Result &res, const char *path) {
    std::ofstream out(path);
    for (int i = 0; i < res.g.size(); ++i) out << res.g(i) << ' ' << res.inv(i) << '\n';
}

}  // namespace

int main() {
    spdlog::set_level(spdlog::level::warn);   // info prints the solver's per-constraint table
    try {
        // The shipped "default" recipe; another one: load_recipe("gropt/diffusion_recipes.json", "fast")
        const Recipe R;

        // ============================================
        // 1) One waveform at a fixed TE
        std::printf("== TE = 80 ms\n");
        const Result fixed = solve_at(80e-3, R);
        report(fixed);
        save(fixed, "diffusion_TE80_10us.txt");

        // ============================================
        // 2) The shortest TE whose waveform reaches b = 1000 s/mm^2
        std::printf("\n== shortest TE for b >= 1000\n");
        const Result shortest = te_search(1000.0, R);
        report(shortest);
        save(shortest, "diffusion_minTE_10us.txt");

        return (fixed.ok && shortest.ok) ? 0 : 1;   // ok: converged and every 10 us check passed
    } catch (const std::exception &e) {             // load_recipe failed, or b not reachable below TE_max
        spdlog::error("{}", e.what());
        return 1;
    }
}
