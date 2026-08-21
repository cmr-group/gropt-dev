// demo_diffusion.cpp -- diffusion solve driven by a RECIPE, the C++ mirror of gropt/diffusion.py.
//
//   gropt                          -> the built-in default recipe (best_pns_1 from the random sweep)
//   gropt recipes.json             -> the first recipe in that library file
//   gropt recipes.json best_pns_1  -> that named recipe
//
// The split matches diffusion_recipes.py: a PROBLEM (timing, hardware limits, which constraints exist --
// what the optimum IS) is hardcoded below, and a RECIPE (weights, x0 seed, solver knobs -- HOW to solve)
// is overlaid on top, either from JSON or from default_recipe(). Every field the Python
// DiffParams/SolverCfg sweep varies is a recipe field here, so a recipes.json written by
// gropt.diffusion_recipes.save_recipe loads with no translation.
#include <cmath>
#include <fstream>
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
#include "gropt_params.hpp"
#include "nlohmann/json.hpp"
#include "op_bvalue.hpp"
#include "solver_groptsdmm.hpp"

using namespace Gropt;

namespace {

#ifdef GROPT_HDF5
// std::vector<Eigen::VectorXd> -> vector<vector<double>> (iters x N) -- a native H5Easy type, so no
// Eigen serialization support is needed. All rows are length N, so it dumps as a rectangular 2D dataset.
std::vector<std::vector<double>> to_rows(const std::vector<Eigen::VectorXd> &v) {
    std::vector<std::vector<double>> out;
    out.reserve(v.size());
    for (const auto &x : v) out.emplace_back(x.data(), x.data() + x.size());
    return out;
}

void save_debug(Solver &solver) {  // needs solver.extra_debug = true to have populated history
    H5Easy::File f("debug_output.h5", H5Easy::File::Overwrite);
    H5Easy::dump(f, "/hist_x", to_rows(solver.debug_solver.hist_X));  // 2D dataset (iters x N)
    H5Easy::dump(f, "/hist_cg_iter", solver.hist_cg_iter);            // std::vector<int>
    spdlog::info("wrote debug_output.h5");
}
#endif

// ===================================================================================================
// SAFE (PNS / cardiac) coefficient tables
// ===================================================================================================
// One SAFE table: 8 vectors of length 3 (axes x,y,z). tau in seconds.
struct SafeCoeffs {
    Eigen::VectorXd tau1{3}, tau2{3}, tau3{3}, a1{3}, a2{3}, a3{3}, stim_limit{3}, g_scale{3};
};

// These are gropt.readasc.get_random_safe_params(42) verbatim -- the SAFE source the Python sweep that
// produced the recipes ran against (SafeSource(kind="random", seed=42), the DiffParams default), so the
// tuned weights here mean what they meant there. Swap in a real scanner table for real work.
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
// PROBLEM -- what the optimum is. Mirrors diffusion_recipes.PROBLEM_FIELDS: never part of a recipe,
// so the same recipe can be reused across geometries/hardware.
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

    // b-value: "obj" maximizes b; "setval"/"minval"/"minval_max" enforce it as a constraint
    std::string bval_mode = "obj";
    double bval_min = 100.0;           // constraint modes only

    // optional constraints (off by default)
    bool concomitant = false;
    double eddy_lam = -1.0;            // [s]; < 0 => off
    double jerk_lam = 0.0;             // order-2 TV weight; 0 => off
    int basin_same_sign = -1;          // -1 off | 0 force sign flip | 1 force same sign
    double basin_window = 1e-3, basin_eps = 0.07;
};

// ===================================================================================================
// RECIPE -- how to solve. Field-for-field the JSON written by diffusion_recipes.save_recipe: the
// "diff" block is the non-problem half of DiffParams, the "solver" block is all of SolverCfg. Defaults
// here match those Python dataclass defaults, so a partial recipe behaves like Python's replace().
// ===================================================================================================
struct Recipe {
    std::string description;

    // --- "diff" block (DiffParams solve knobs) ---
    double w_gmax = 1.0, w_smax = 1.0, w_moment = 1.0;
    double w_pns = 1.0, w_cns = 1.0, w_concomitant = 1.0, w_eddy = 1.0, w_jerk = 1.0, w_bval = 1.0;
    bool   moment_project = true;
    double safe_eps = 0.0;
    double bval_obj_weight = 1.0;   // obj mode
    double bval_max_scale = 1.02;   // constraint mode
    std::string x0_mode = "diff_init";   // "diff_init" | "const" | "sine"
    double x0_amp = 0.01;
    bool   x0_invert = true;
    double x0_periods = 1.0;
    bool   x0_project = false;

    // --- "solver" block (SolverCfg) ---
    int    max_iter = 4000, max_feval = 200000, min_iter = 1, obj_patience = 20;
    double obj_rtol = 1e-4, gamma_x = 1.6;

    double ils_tol = 0.1;
    int    ils_max_iter = 20, ils_min_iter = 2;
    double ils_sigma = 1e-4, ils_tik_lam = 0.0;

    bool   bb_reweight = true;      // BB per-operator adaptation
    int    rw_interval = 8;
    double rw_e_corr = 0.2, rw_scalelim = 2.0, rw_eps = 1e-36;

    bool   grw = false;             // global "double the most-infeasible op"
    int    grw_interval = 20;
    double grw_mod = 2.0;
    bool   grw_balanced = false;

    bool   reproject_iterate = true;

    double cutoff_freq = -1.0;
    int    cutoff_iter = -1;
    double cutoff_trans = 0.0;

    bool   tr_enable = false;       // trust-region step control
    double tr_tol = -1.0, tr_bump = 4.0;
    int    tr_max_reject = 5;
    double tr_decay = 0.5;
    std::string tr_monitor = "linearization_error";

    bool   obj_gate = false;        // feasibility-gated objective
    double obj_gate_scale = 0.05;

    bool   extra_debug = false;
};

// The tuned sweep winner (recipes.json "best_pns_1", frac_optimal=0.983), inline so the demo runs
// with the good settings without a JSON file.
Recipe default_recipe() {
    Recipe R;
    R.description = "built-in best_pns_1 (frac_optimal=0.983)";

    R.w_smax = 312.5601580523615;
    R.w_pns = 47.05829778348828;
    R.w_cns = 10.547056146456551;
    R.bval_obj_weight = 0.21721639516550184;
    R.bval_max_scale = 1.02;
    R.x0_mode = "sine";
    R.x0_amp = 0.036865715905310396;
    R.x0_invert = false;
    R.x0_periods = 4;
    R.x0_project = true;

    R.max_iter = 4000;
    R.gamma_x = 1.0816673248634217;
    R.ils_tol = 0.030682565018240758;
    R.ils_max_iter = 26;
    R.ils_min_iter = 2;
    R.ils_sigma = 1.5214267699983488e-05;
    R.rw_interval = 33;
    R.rw_e_corr = 0.04925816096961088;
    R.rw_scalelim = 9.591148924091835;
    R.grw = true;
    R.grw_interval = 15;
    R.grw_mod = 9.56506967232672;
    R.grw_balanced = false;
    R.reproject_iterate = true;
    R.tr_enable = true;
    R.tr_bump = 10.985469703019819;
    R.tr_tol = 0.1356731687978308;
    R.tr_decay = 0.22704167936390063;
    R.tr_max_reject = 25;
    R.tr_monitor = "linearization_error";
    R.obj_gate = false;
    R.obj_gate_scale = 0.05;
    return R;
}

// ===================================================================================================
// Recipe JSON
// ===================================================================================================
// Overlay one JSON block onto R. Unknown keys warn instead of failing, so a recipe written by a newer
// Python (a knob this build has no setter for) still loads -- but you are told what was dropped.
void apply_diff_block(const nlohmann::json &d, Recipe &R) {
    for (const auto &[k, v] : d.items()) {
        try {
            if      (k == "w_gmax")          R.w_gmax = v.get<double>();
            else if (k == "w_smax")          R.w_smax = v.get<double>();
            else if (k == "w_moment")        R.w_moment = v.get<double>();
            else if (k == "w_pns")           R.w_pns = v.get<double>();
            else if (k == "w_cns")           R.w_cns = v.get<double>();
            else if (k == "w_concomitant")   R.w_concomitant = v.get<double>();
            else if (k == "w_eddy")          R.w_eddy = v.get<double>();
            else if (k == "w_jerk")          R.w_jerk = v.get<double>();
            else if (k == "w_bval")          R.w_bval = v.get<double>();
            else if (k == "moment_project")  R.moment_project = v.get<bool>();
            else if (k == "safe_eps")        R.safe_eps = v.get<double>();
            else if (k == "bval_obj_weight") R.bval_obj_weight = v.get<double>();
            else if (k == "bval_max_scale")  R.bval_max_scale = v.get<double>();
            else if (k == "x0_mode")         R.x0_mode = v.get<std::string>();
            else if (k == "x0_amp")          R.x0_amp = v.get<double>();
            else if (k == "x0_invert")       R.x0_invert = v.get<bool>();
            else if (k == "x0_periods")      R.x0_periods = v.get<double>();
            else if (k == "x0_project")      R.x0_project = v.get<bool>();
            else spdlog::warn("recipe: ignoring unknown 'diff' key '{}'", k);
        } catch (const nlohmann::json::exception &e) {   // nlohmann's type error does not name the key
            throw std::runtime_error("recipe 'diff' key '" + k + "': " + e.what());
        }
    }
}

void apply_solver_block(const nlohmann::json &d, Recipe &R) {
    for (const auto &[k, v] : d.items()) {
        try {
            if      (k == "max_iter")           R.max_iter = v.get<int>();
            else if (k == "max_feval")          R.max_feval = v.get<int>();
            else if (k == "min_iter")           R.min_iter = v.get<int>();
            else if (k == "obj_patience")       R.obj_patience = v.get<int>();
            else if (k == "obj_rtol")           R.obj_rtol = v.get<double>();
            else if (k == "gamma_x")            R.gamma_x = v.get<double>();
            else if (k == "ils_tol")            R.ils_tol = v.get<double>();
            else if (k == "ils_max_iter")       R.ils_max_iter = v.get<int>();
            else if (k == "ils_min_iter")       R.ils_min_iter = v.get<int>();
            else if (k == "ils_sigma")          R.ils_sigma = v.get<double>();
            else if (k == "ils_tik_lam")        R.ils_tik_lam = v.get<double>();
            else if (k == "bb_reweight")        R.bb_reweight = v.get<bool>();
            else if (k == "rw_interval")        R.rw_interval = v.get<int>();
            else if (k == "rw_e_corr")          R.rw_e_corr = v.get<double>();
            else if (k == "rw_scalelim")        R.rw_scalelim = v.get<double>();
            else if (k == "rw_eps")             R.rw_eps = v.get<double>();
            else if (k == "grw")                R.grw = v.get<bool>();
            else if (k == "grw_interval")       R.grw_interval = v.get<int>();
            else if (k == "grw_mod")            R.grw_mod = v.get<double>();
            else if (k == "grw_balanced")       R.grw_balanced = v.get<bool>();
            else if (k == "reproject_iterate")  R.reproject_iterate = v.get<bool>();
            else if (k == "cutoff_freq")        R.cutoff_freq = v.get<double>();
            else if (k == "cutoff_iter")        R.cutoff_iter = v.get<int>();
            else if (k == "cutoff_trans")       R.cutoff_trans = v.get<double>();
            else if (k == "tr_enable")          R.tr_enable = v.get<bool>();
            else if (k == "tr_tol")             R.tr_tol = v.get<double>();
            else if (k == "tr_bump")            R.tr_bump = v.get<double>();
            else if (k == "tr_max_reject")      R.tr_max_reject = v.get<int>();
            else if (k == "tr_decay")           R.tr_decay = v.get<double>();
            else if (k == "tr_monitor")         R.tr_monitor = v.get<std::string>();
            else if (k == "obj_gate")           R.obj_gate = v.get<bool>();
            else if (k == "obj_gate_scale")     R.obj_gate_scale = v.get<double>();
            else if (k == "extra_debug")        R.extra_debug = v.get<bool>();
            else spdlog::warn("recipe: ignoring unknown 'solver' key '{}'", k);
        } catch (const nlohmann::json::exception &e) {   // nlohmann's type error does not name the key
            throw std::runtime_error("recipe 'solver' key '" + k + "': " + e.what());
        }
    }
}

// Load `name` (empty => the first entry) from a recipe library JSON. Missing keys keep the Recipe
// default, exactly like Python's replace(base_scfg, **r["solver"]).
// nlohmann::json keys are sorted, which is also how save_recipe writes them (sort_keys=True), so
// "the first entry" is both the alphabetically first and the first in a Python-written file.
Recipe load_recipe(const std::string &path, const std::string &name) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("could not open recipe file: " + path);

    nlohmann::json lib;
    try {
        f >> lib;
    } catch (const nlohmann::json::parse_error &e) {
        throw std::runtime_error("could not parse " + path + ": " + e.what());
    }
    if (!lib.is_object() || lib.empty()) {
        throw std::runtime_error("recipe library is not a non-empty JSON object: " + path);
    }

    auto it = name.empty() ? lib.begin() : lib.find(name);
    if (it == lib.end()) {
        std::string have;
        for (const auto &el : lib.items()) have += (have.empty() ? "" : ", ") + el.key();
        throw std::runtime_error("recipe '" + name + "' not found in " + path + " (have: " + have + ")");
    }
    const nlohmann::json &entry = it.value();

    Recipe R;
    if (entry.contains("description")) R.description = entry["description"].get<std::string>();
    if (entry.contains("diff")) apply_diff_block(entry["diff"], R);
    if (entry.contains("solver")) apply_solver_block(entry["solver"], R);
    spdlog::info("loaded recipe '{}' from {} ({})", it.key(), path, R.description);
    return R;
}

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

// Build the problem into `gp` (an out-param: GroptParams holds references into its own pdata, so it
// must not be moved/returned by value). Returns start_idx (> 0 only for preencode).
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
        if (P.pns_lim >= 0.0) {
            const SafeCoeffs s = pns_table();
            gp.add_SAFE(P.pns_lim, s.tau1, s.tau2, s.tau3, s.a1, s.a2, s.a3, s.stim_limit, s.g_scale, 0, R.w_pns);
        }
        if (P.cns_lim >= 0.0) {
            const SafeCoeffs s = cns_table();
            gp.add_SAFE(P.cns_lim, s.tau1, s.tau2, s.tau3, s.a1, s.a2, s.a3, s.stim_limit, s.g_scale, 0, R.w_cns);
        }
    }

    if (P.concomitant) gp.add_concomitant(start_idx, true, R.w_concomitant, 0.1, 1.0, false, 1.0, true);
    if (P.eddy_lam >= 0.0) {
        Eigen::VectorXd lam(1);
        lam(0) = P.eddy_lam;
        gp.add_eddy(lam, 1e-4, R.w_eddy, true);
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

    solver.grw_balanced = R.grw_balanced;

    solver.set_ils_params(R.ils_tol, R.ils_max_iter, R.ils_min_iter, R.ils_sigma, R.ils_tik_lam);

    // Reweighter gating, same trick as diffusion.py make_solver (there is no do_rw setter):
    //   BB is gated by `iiter > rw_interval` -> push rw_interval past max_iter to disable it.
    //   grw multiplies by grw_mod -> grw_mod = 1.0 is a no-op.
    const int rw_interval = R.bb_reweight ? R.rw_interval : (R.max_iter + 1);
    const double grw_mod = R.grw ? R.grw_mod : 1.0;
    solver.set_sdmm_params(rw_interval, R.rw_e_corr, R.rw_eps, R.rw_scalelim, 20, R.grw_interval, grw_mod);
}

}  // namespace

void demo_diffusion(const std::string &recipe_path, const std::string &recipe_name) {
    const Problem P;
    const Recipe R = recipe_path.empty() ? default_recipe() : load_recipe(recipe_path, recipe_name);

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
