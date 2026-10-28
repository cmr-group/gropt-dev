// demo_diffusion.cpp -- recipe-driven diffusion solve; the C++ counterpart of gropt/diffusion.py.
//
//   gropt                                   -> the built-in default recipe (Recipe{})
//   gropt gropt/diffusion_recipes.json fast -> a role or entry of the shipped library
//   gropt my_recipes.json [name]            -> a save_recipe file (default: its first recipe)
//
// The Problem (timing, limits, constraints) is fixed below; the Recipe (weights, x0 seed, solver settings)
// comes from diffusion_recipe.hpp.
#include <stdexcept>
#include <string>
#include <vector>

#include "spdlog/spdlog.h"
#ifdef GROPT_PLOTTING
#include <matplot/matplot.h>
#endif
#ifdef GROPT_HDF5
#include <highfive/H5Easy.hpp>
#include "solver.hpp"
#endif
#include "diffusion_recipe.hpp"
#include "gropt_params.hpp"
#include "op_bvalue.hpp"
#include "op_safe_slack.hpp"
#include "solver_groptsdmm.hpp"

using namespace Gropt;

namespace {

#ifdef GROPT_HDF5
// Equal-length rows, so H5Easy writes an (iters x N) 2D dataset without Eigen support.
std::vector<std::vector<double>> to_rows(const std::vector<Eigen::VectorXd> &v) {
    std::vector<std::vector<double>> out;
    out.reserve(v.size());
    for (const auto &x : v) out.emplace_back(x.data(), x.data() + x.size());
    return out;
}

void save_debug(Solver &solver) {  // needs solver.extra_debug = true to have populated history
    H5Easy::File f("debug_output.h5", H5Easy::File::Overwrite);
    H5Easy::dump(f, "/hist_x", to_rows(solver.debug_solver.hist_X));
    H5Easy::dump(f, "/hist_cg_iter", solver.debug_solver.hist_cg_iter);
    spdlog::info("wrote debug_output.h5");
}
#endif

// ===================================================================================================
// SAFE (PNS / cardiac) coefficient tables
// ===================================================================================================
// One entry per axis (x, y, z); tau [s].
struct SafeCoeffs {
    Eigen::VectorXd tau1{3}, tau2{3}, tau3{3}, a1{3}, a2{3}, a3{3}, stim_limit{3}, g_scale{3};
};

// gropt.readasc.get_random_safe_params(42), the DiffParams default, so Python-tuned recipes carry over.
// Use a real scanner table (e.g. from an .asc file) for real work.
SafeCoeffs pns_table() {
    SafeCoeffs s;
    s.tau1 << 0.86e-3, 0.93e-3, 0.78e-3;
    s.tau2 << 11.0e-3, 13.0e-3, 10.0e-3;
    s.tau3 << 0.27e-3, 0.23e-3, 0.25e-3;
    s.a1   << 0.28, 0.24, 0.29;
    s.a2   << 0.54, 0.42, 0.60;
    s.a3   << 0.18, 0.34, 0.11;
    s.stim_limit << 29.0, 27.4, 38.5;
    s.g_scale    << 0.34, 0.34, 0.29;
    return s;
}

SafeCoeffs cns_table() {
    SafeCoeffs s;
    s.tau1 << 2.7e-3, 3.0e-3, 2.3e-3;
    s.tau2 << 1.4e-3, 1.5e-3, 1.2e-3;
    s.tau3 << 1.1e-3, 1.5e-3, 1.2e-3;
    s.a1   << 0.67, 0.79, 0.78;
    s.a2   << 0.33, 0.21, 0.22;
    s.a3   << 0.00, 0.00, 0.00;
    s.stim_limit << 14.3, 14.9, 18.1;
    s.g_scale    << 0.34, 0.30, 0.32;
    return s;
}

// ===================================================================================================
// Problem: timing, hardware limits and constraints (diffusion_recipes.PROBLEM_FIELDS; never in a recipe)
// ===================================================================================================
struct Problem {
    // timing [s]
    std::string diff_mode = "gropt";   // "gropt" | "conventional" | "preencode"
    double dt = 400e-6, TE = 80e-3, T_90 = 3e-3, T_180 = 5e-3, T_readout = 16e-3;
    double T_pre = 0.0;                // preencode only

    // hardware limits
    double gmax = 0.19;                // [T/m]
    double smax = 200.0;               // [T/m/s]

    // moments: null M0..M_MMT
    int MMT = 1;
    double moment_tol = 1e-4;

    // SAFE thresholds; < 0 => that model is off
    double pns_lim = 0.8;
    double cns_lim = 0.8;
    bool safe_alpha_exact = true;      // raster-invariant filter (what the recipes are tuned with)

    // b-value: "obj" maximizes b; "setval"/"minval"/"minval_max" enforce it as a constraint
    std::string bval_mode = "obj";
    double bval_min = 100.0;           // constraint modes only

    // optional constraints (off by default)
    bool concomitant = false;
    double concomitant_tol = 0.1;      // max |E_pre/E_post - 1|
    bool concomitant_project = true;   // exact balance; false = the soft band
    double eddy_lam = -1.0;            // [s]; < 0 => off
    double jerk_lam = 0.0;             // order-2 TV weight; 0 => off
    int basin_same_sign = -1;          // -1 off | 0 force sign flip | 1 force same sign
    double basin_window = 1e-3, basin_eps = 0.07;
};

// ===================================================================================================
// Builder (the C++ mirror of diffusion.py's build_gparams; the seed and solver settings are in
// diffusion_recipe.hpp)
// ===================================================================================================
// Build into `gp` (out-param: GroptParams holds references into its own pdata, so it can't be returned
// by value). Returns start_idx (> 0 only for preencode).
int build_gparams(const Problem &P, const Recipe &R, GroptParams &gp) {
    int start_idx = 0;
    if (P.diff_mode == "gropt") {
        gp.diff_init(P.dt, P.TE, P.T_90, P.T_180, P.T_readout);
    } else if (P.diff_mode == "conventional") {
        gp.diff_init_deadtime(P.dt, P.TE, P.T_90, P.T_180, P.T_readout);
    } else if (P.diff_mode == "preencode") {
        start_idx = gp.diff_init_preencode(P.dt, P.TE, P.T_90, P.T_180, P.T_readout, P.T_pre);
    } else {
        throw std::invalid_argument("Unknown diff_mode: " + P.diff_mode);
    }

    gp.add_gmax(P.gmax, true, R.w_gmax);
    gp.add_smax(P.smax, true, R.w_smax);
    for (int m = 0; m <= P.MMT; ++m)                    // null M0..M_MMT
        gp.add_moment(m, 0.0, P.moment_tol, "mT*ms/m", 0, start_idx, -1, 0, R.w_moment, R.moment_project);

    if (P.pns_lim >= 0.0 || P.cns_lim >= 0.0) {
        gp.safe_eps = R.safe_eps;                       // copied into each Op_SAFE by add_SAFE; set first
        gp.safe_signed13 = R.safe_signed13;
        gp.safe_lifted = R.safe_lifted;
        gp.safe_alpha_exact = P.safe_alpha_exact;
        if (R.safe_lifted) gp.ensure_slack_abs(SAFE_SLACK_BLOCK, R.w_slack);   // first, as diffusion.py does
        if (P.pns_lim >= 0.0) {
            const SafeCoeffs s = pns_table();
            gp.add_SAFE(P.pns_lim, s.tau1, s.tau2, s.tau3, s.a1, s.a2, s.a3, s.stim_limit, s.g_scale, 0, R.w_pns);
        }
        if (P.cns_lim >= 0.0) {
            const SafeCoeffs s = cns_table();
            gp.add_SAFE(P.cns_lim, s.tau1, s.tau2, s.tau3, s.a1, s.a2, s.a3, s.stim_limit, s.g_scale, 0, R.w_cns);
        }
    }

    if (P.concomitant)
        gp.add_concomitant(start_idx, true, R.w_concomitant, P.concomitant_tol, 1.0, P.concomitant_project);
    if (P.eddy_lam >= 0.0) {
        Eigen::VectorXd lam(1);
        lam(0) = P.eddy_lam;
        gp.add_eddy(lam, 1e-4, R.w_eddy, R.eddy_project);
    }
    if (P.jerk_lam > 0.0) gp.add_TV(P.jerk_lam, R.w_jerk, 2);
    if (P.basin_same_sign >= 0)
        gp.add_diff_basin(P.basin_window, P.basin_eps, P.gmax, 1.0, P.basin_same_sign == 1);

    // b-value: maximize (objective) or enforce (constraint). The objective path ignores
    // target/tol/mode/max_scale; these mirror the Python add_bvalue defaults.
    if (P.bval_mode == "obj") {
        gp.add_bvalue(100.0, 1.0, start_idx, -1, R.bval_obj_weight, BVALUE_MODE_MINVAL, 1.01, true, true);
        gp.normalize_obj = true;
    } else {
        const int mode = (P.bval_mode == "setval")     ? BVALUE_MODE_SETVAL
                       : (P.bval_mode == "minval")     ? BVALUE_MODE_MINVAL
                       : (P.bval_mode == "minval_max") ? BVALUE_MODE_MINVALMAX
                                                       : -1;
        if (mode < 0) throw std::invalid_argument("Unknown bval_mode: " + P.bval_mode);
        gp.add_bvalue(P.bval_min, 1.0, start_idx, -1, R.w_bval, mode, R.bval_max_scale, false, true);
    }

    set_recipe_x0(gp, R, P.MMT, start_idx);

    gp.prepare();
    return start_idx;
}

}  // namespace

void demo_diffusion(const std::string &recipe_path, const std::string &recipe_name) {
    const Problem P;
    const Recipe R = recipe_path.empty() ? Recipe{} : load_recipe(recipe_path, recipe_name);

    GroptParams gp;
    build_gparams(P, R, gp);

    SolverGroptSDMM solver;
    configure_solver(R, solver);

#ifdef GROPT_HDF5
    solver.extra_debug = true;   // capture per-iteration history for save_debug
#endif

    SolveResult res = solver.solve(gp);
    spdlog::info("b = {:.1f}  converged = {}  iters = {}", res.bvalue, res.converged, res.n_iter);

#ifdef GROPT_HDF5
    save_debug(solver);
#endif

#ifdef GROPT_PLOTTING
    std::vector<double> g(res.X.data(), res.X.data() + res.X.size());
    matplot::plot(g);
    matplot::title("diffusion gradient (PNS + CNS)");
    matplot::xlabel("sample");
    matplot::ylabel("G [T/m]");
    matplot::show();
#endif
}
