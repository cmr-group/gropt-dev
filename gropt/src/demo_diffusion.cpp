// demo_diffusion.cpp -- recipe-driven diffusion solve; the C++ counterpart of gropt/diffusion.py.
//
//   gropt                                   -> the built-in default recipe (Recipe{})
//   gropt gropt/diffusion_recipes.json fast -> a role or entry of the shipped library
//   gropt my_recipes.json [name]            -> a save_recipe file (default: its first recipe)
//
// The Problem (timing, limits, constraints) is fixed below; the Recipe (weights, x0 seed, solver settings)
// comes from diffusion_recipe.hpp.
#include <cmath>
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
// Builders (the C++ mirror of diffusion.py's _build_x0 / build_gparams / make_solver)
// ===================================================================================================
// X0 seed. Empty return -> keep the diff_init seed. Call after diff_init, before prepare.
Eigen::VectorXd build_x0(const GroptParams &gp, const Recipe &R, int MMT, int start_idx) {
    if (R.x0_mode == "diff_init") return {};
    const int N = gp.N;
    const double dt = gp.dt;
    const Eigen::VectorXd &inv = gp.pdata.inv_vec, &fixer = gp.pdata.fixer, &setv = gp.pdata.set_vals;
    auto free = [&](int i) { return fixer(i) > 0.5; };  // 1 = free, 0 = fixed
    constexpr double PI = 3.14159265358979323846;
    Eigen::VectorXd x = Eigen::VectorXd::Zero(N);

    if (R.x0_mode == "const") {
        for (int i = 0; i < N; ++i) if (free(i)) x(i) = R.x0_amp;
        if (R.x0_invert) x.array() *= inv.array();
    } else if (R.x0_mode == "sine") {
        for (int i = 0; i < N;) {                       // full sine per free run: 0 at both ends
            if (!free(i)) { ++i; continue; }
            int j = i; while (j < N && free(j)) ++j;
            const int n = j - i;
            for (int k = 0; k < n && n >= 2; ++k)
                x(i + k) = R.x0_amp * std::sin(2.0 * PI * R.x0_periods * (double(k) / (n - 1)));
            i = j;
        }
        if (R.x0_invert) for (int i = 0; i < N; ++i) if (inv(i) < 0) x(i) = -x(i);
    } else {
        throw std::invalid_argument("Unknown x0_mode: " + R.x0_mode);
    }

    for (int i = 0; i < N; ++i) if (!free(i)) x(i) = std::isnan(setv(i)) ? 0.0 : setv(i);

    if (R.x0_project && MMT >= 0) {                      // project onto moment null-space (free DOFs)
        const int nm = MMT + 1;
        Eigen::MatrixXd M(nm, N), Mt(nm, N);
        for (int k = 0; k < nm; ++k)
            for (int i = 0; i < N; ++i) {
                M(k, i)  = (i < start_idx) ? 0.0 : std::pow(i * dt, k) * inv(i);
                Mt(k, i) = M(k, i) * fixer(i);
            }
        x -= Mt.transpose() * (Mt * Mt.transpose()).ldlt().solve(M * x);
    }
    return x;
}

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

    Eigen::VectorXd x0 = build_x0(gp, R, P.MMT, start_idx);   // empty -> keep the diff_init seed
    if (x0.size()) gp.setvec_X0(x0, 1, false);

    gp.prepare();
    return start_idx;
}

void configure_solver(const Recipe &R, SolverGroptSDMM &solver) {
    solver.max_iter = R.max_iter;
    solver.max_feval = R.max_feval;
    solver.min_iter = R.min_iter;
    solver.obj_patience = R.obj_patience;
    solver.obj_rtol = R.obj_rtol;
    solver.gamma_x = R.gamma_x;
    solver.extra_debug = R.extra_debug;

    solver.reproject_iterate = R.reproject_iterate;

    solver.cutoff_freq = R.cutoff_freq;
    solver.cutoff_iter = R.cutoff_iter;
    solver.cutoff_trans = R.cutoff_trans;

    solver.tr_enable = R.tr_enable;
    solver.tr_tol = R.tr_tol;
    solver.tr_bump = R.tr_bump;
    solver.tr_max_reject = R.tr_max_reject;
    solver.tr_decay = R.tr_decay;
    solver.tr_monitor = R.tr_monitor;

    solver.obj_gate_enable = R.obj_gate;
    solver.obj_gate_scale = R.obj_gate_scale;

    solver.ils_tol = R.ils_tol;
    solver.ils_max_iter = R.ils_max_iter;
    solver.ils_min_iter = R.ils_min_iter;
    solver.ils_sigma = R.ils_sigma;
    solver.ils_tik_lam = R.ils_tik_lam;

    solver.bb_enable = R.bb_reweight;
    solver.rw_interval = R.rw_interval;
    solver.rw_e_corr = R.rw_e_corr;
    solver.rw_eps = R.rw_eps;
    solver.rw_scalelim = R.rw_scalelim;
    solver.grw_enable = R.grw;
    solver.grw_interval = R.grw_interval;
    solver.grw_mod = R.grw_mod;
    solver.grw_balanced = R.grw_balanced;
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
