#ifndef DIFFUSION_RECIPE_H
#define DIFFUSION_RECIPE_H

// Diffusion solve recipe, the C++ mirror of gropt/diffusion_recipes.py: the solve knobs only ("diff" =
// DiffParams' RECIPE_FIELDS, "solver" = SolverCfg); the problem stays with the caller. The defaults are the
// tuned "default" recipe of gropt/diffusion_recipes.json (tests/test_diffusion_recipe.cpp checks they match).

#include <fstream>
#include <set>
#include <stdexcept>
#include <string>

#include "nlohmann/json.hpp"
#include "spdlog/spdlog.h"

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

}  // namespace Gropt

#endif
