#ifndef DIFFUSION_RECIPE_H
#define DIFFUSION_RECIPE_H

// Diffusion solve recipe, the C++ mirror of gropt/diffusion_recipes.py: the solve knobs only ("diff" =
// DiffParams' RECIPE_FIELDS, "solver" = SolverCfg); the problem stays with the caller. The defaults are the
// tuned "default" recipe of gropt/diffusion_recipes.json (tests/test_diffusion_recipe.cpp checks they match).
//
// Using one: pass the R.w_* weights to the add_* calls and copy R.safe_* to the GroptParams before add_SAFE,
// then set_recipe_x0() before prepare() and configure_solver() on the solver (see demo_diffusion_simple.cpp).

#include <cmath>
#include <fstream>
#include <set>
#include <stdexcept>
#include <string>

#include "nlohmann/json.hpp"
#include "spdlog/spdlog.h"

#include "gropt_params.hpp"
#include "solver_groptsdmm.hpp"

namespace Gropt {

// X(type, name, default), one per field, so the struct, the reader and the writer cannot drift apart.
#define GROPT_RECIPE_DIFF(X)                                                                                    \
    X(double, w_gmax, 1.0) X(double, w_smax, 10.832717264505106) X(double, w_moment, 1.0)                       \
    X(double, w_bval, 1.0) X(double, w_pns, 12.366370987704867) X(double, w_cns, 57.339230929159065)            \
    X(double, w_slack, 105.08823640819529) X(double, w_concomitant, 4.371489218804275) X(double, w_eddy, 1.0)   \
    X(double, w_jerk, 1.0) X(bool, moment_project, true) X(bool, eddy_project, true)                            \
    X(double, bval_obj_weight, 1.7034494007148022) X(double, bval_max_scale, 1.02)                              \
    X(double, safe_eps, 0.0) X(bool, safe_signed13, false) X(bool, safe_lifted, true)                           \
    X(std::string, x0_mode, "sine") X(double, x0_amp, 0.08443255496637647) X(bool, x0_invert, false)            \
    X(double, x0_periods, 2.0) X(bool, x0_project, false)

#define GROPT_RECIPE_SOLVER(X)                                                                                  \
    X(int, max_iter, 4000) X(int, max_feval, 200000) X(int, min_iter, 1) X(int, obj_patience, 20)               \
    X(double, obj_rtol, 1e-4) X(double, gamma_x, 1.165149741846207)                                             \
    X(double, ils_tol, 0.07290432029478086) X(int, ils_max_iter, 31) X(int, ils_min_iter, 2)                    \
    X(double, ils_sigma, 0.0001554544663042038) X(double, ils_tik_lam, 0.0)                                     \
    X(bool, bb_reweight, true) X(int, rw_interval, 37) X(double, rw_e_corr, 0.16079028465931097)                \
    X(double, rw_scalelim, 4.991259555225381) X(double, rw_eps, 1e-36)                                          \
    X(bool, grw, true) X(int, grw_interval, 31) X(double, grw_mod, 6.005900533768388) X(bool, grw_balanced, false) \
    X(bool, reproject_iterate, false)                                                                           \
    X(double, cutoff_freq, -1.0) X(int, cutoff_iter, -1) X(double, cutoff_trans, 0.0)                           \
    X(bool, tr_enable, false) X(double, tr_tol, -1.0) X(double, tr_bump, 4.0) X(int, tr_max_reject, 5)         \
    X(double, tr_decay, 0.5) X(std::string, tr_monitor, "linearization_error")                                  \
    X(bool, obj_gate, false) X(double, obj_gate_scale, 0.05) X(bool, extra_debug, false)

struct Recipe {
    std::string description = "built-in default";
#define GROPT_RECIPE_FIELD(type, name, dflt) type name = dflt;
    GROPT_RECIPE_DIFF(GROPT_RECIPE_FIELD)
    GROPT_RECIPE_SOLVER(GROPT_RECIPE_FIELD)
#undef GROPT_RECIPE_FIELD

    // {"diff": ..., "solver": ...}, the blocks save_recipe writes
    nlohmann::json to_json() const {
        nlohmann::json d, s;
#define GROPT_RECIPE_DIFF_OUT(type, name, dflt) d[#name] = name;
#define GROPT_RECIPE_SOLVER_OUT(type, name, dflt) s[#name] = name;
        GROPT_RECIPE_DIFF(GROPT_RECIPE_DIFF_OUT)
        GROPT_RECIPE_SOLVER(GROPT_RECIPE_SOLVER_OUT)
#undef GROPT_RECIPE_DIFF_OUT
#undef GROPT_RECIPE_SOLVER_OUT
        return {{"diff", d}, {"solver", s}};
    }
};

// Overlay one "diff" or "solver" block onto R. An unknown key is an error: a knob this build cannot apply
// would silently change the solve. Older files stored some fields that are now problem fields; skip those.
inline void apply_recipe_block(const nlohmann::json &blk, const std::string &which, Recipe &R) {
    static const std::set<std::string> problem_keys{"concomitant_project", "safe_alpha_exact",
                                                    "concomitant_exact_quad", "moment_pwl_quad", "bval_pwl_quad"};
    for (const auto &[k, v] : blk.items()) {
        try {
#define GROPT_RECIPE_READ(type, name, dflt) if (k == #name) { R.name = v.get<type>(); continue; }
            if (which == "diff") { GROPT_RECIPE_DIFF(GROPT_RECIPE_READ) }
            else { GROPT_RECIPE_SOLVER(GROPT_RECIPE_READ) }
#undef GROPT_RECIPE_READ
        } catch (const nlohmann::json::exception &e) {   // nlohmann's type error does not name the key
            throw std::runtime_error("recipe '" + which + "' key '" + k + "': " + e.what());
        }
        if (which == "diff" && problem_keys.count(k)) continue;
        throw std::runtime_error("recipe '" + which + "' has unknown key '" + k + "' (saved by a newer gropt?)");
    }
}

// Load a recipe from a library file: the shipped gropt/diffusion_recipes.json (a role such as "default" or
// "fast", or a dated entry) or a save_recipe file. Empty name: "default", or a save_recipe file's first entry.
// Portfolios (several solves, best verified result) run from Python only.
inline Recipe load_recipe(const std::string &path, const std::string &name) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("could not open recipe file: " + path);
    nlohmann::json lib;
    try {
        f >> lib;
    } catch (const nlohmann::json::parse_error &e) {
        throw std::runtime_error("could not parse " + path + ": " + e.what());
    }
    if (!lib.is_object() || lib.empty()) throw std::runtime_error("not a recipe library: " + path);

    const bool curated = lib.value("schema", 0) == 2;
    const nlohmann::json &recipes = curated ? lib.at("recipes") : lib;
    std::string key = name;
    if (curated) {
        if (lib.contains("portfolios") && lib["portfolios"].contains(name))
            throw std::runtime_error("'" + name + "' is a portfolio; run it from Python (gropt.diffusion.solve)");
        if (key.empty()) key = "default";
        if (lib.contains("roles") && lib["roles"].contains(key)) key = lib["roles"][key].get<std::string>();
    } else if (key.empty()) {
        key = recipes.begin().key();
    }
    if (!recipes.contains(key)) {
        std::string have;
        for (const auto &el : (curated ? lib.at("roles") : recipes).items())
            have += (have.empty() ? "" : ", ") + el.key();
        throw std::runtime_error("recipe '" + name + "' not found in " + path + " (have: " + have + ")");
    }
    const nlohmann::json &entry = recipes.at(key);

    Recipe R;
    R.safe_lifted = false;   // a recipe saved before the lifting existed was tuned unlifted
    R.description = entry.value("description", key);
    if (entry.contains("diff")) apply_recipe_block(entry["diff"], "diff", R);
    if (entry.contains("solver")) apply_recipe_block(entry["solver"], "solver", R);
    spdlog::info("loaded recipe '{}' ({}) from {}", name.empty() ? key : name, key, path);
    return R;
}

// The recipe's seed waveform (diffusion.py's _build_x0). Call after diff_init and the add_* calls, before
// prepare(). x0_mode "diff_init" keeps diff_init's seed. MMT: highest nulled moment, for x0_project.
inline void set_recipe_x0(GroptParams &gp, const Recipe &R, int MMT, int start_idx = 0) {
    if (R.x0_mode == "diff_init") return;
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
    gp.setvec_X0(x, 1, false);
}

// The recipe's solver settings (diffusion.py's make_solver).
inline void configure_solver(const Recipe &R, SolverGroptSDMM &solver) {
    // Iteration limits; stop once b gains < obj_rtol over obj_patience feasible iterations
    solver.max_iter = R.max_iter;
    solver.max_feval = R.max_feval;
    solver.min_iter = R.min_iter;
    solver.obj_patience = R.obj_patience;
    solver.obj_rtol = R.obj_rtol;
    solver.gamma_x = R.gamma_x;            // over-relaxation of the x update
    solver.extra_debug = R.extra_debug;

    solver.reproject_iterate = R.reproject_iterate;   // re-project onto the exact equalities every iteration

    // Optional low-pass of the iterate (off in the shipped recipes)
    solver.cutoff_freq = R.cutoff_freq;
    solver.cutoff_iter = R.cutoff_iter;
    solver.cutoff_trans = R.cutoff_trans;

    // Optional trust region on each step (off in the shipped recipes)
    solver.tr_enable = R.tr_enable;
    solver.tr_tol = R.tr_tol;
    solver.tr_bump = R.tr_bump;
    solver.tr_max_reject = R.tr_max_reject;
    solver.tr_decay = R.tr_decay;
    solver.tr_monitor = R.tr_monitor;

    // Optional: damp the objective while SAFE is violated (off in the shipped recipes)
    solver.obj_gate_enable = R.obj_gate;
    solver.obj_gate_scale = R.obj_gate_scale;

    // Inner linear solve (conjugate gradient) per iteration
    solver.ils_tol = R.ils_tol;
    solver.ils_max_iter = R.ils_max_iter;
    solver.ils_min_iter = R.ils_min_iter;
    solver.ils_sigma = R.ils_sigma;
    solver.ils_tik_lam = R.ils_tik_lam;

    // Constraint weight adaptation: per-constraint (bb_*/rw_*) and global (grw_*)
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

}  // namespace Gropt

#endif
