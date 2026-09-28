// The C++ recipe reader (diffusion_recipe.hpp) against the shipped library: Recipe{} is its "default", every
// role loads, portfolios and unknown keys are refused, and problem keys in older files are skipped.

#include "diffusion_recipe.hpp"
#include "test_util.hpp"

#include <cstdio>
#include <filesystem>
#include <fstream>

using namespace Gropt;
using namespace gropt_test;

namespace {

template <class F> bool throws_runtime(F &&f) {
    try {
        f();
    } catch (const std::runtime_error &) {
        return true;
    }
    return false;
}

}  // namespace

int run_diffusion_recipe_tests() {
    std::printf("\nDiffusion recipes (C++ reader)\n");
    int failures = 0;
    const std::string lib_path = GROPT_RECIPE_LIBRARY;

    failures += report(load_recipe(lib_path, "default").to_json() == Recipe{}.to_json(),
                       "Recipe{} == the shipped library's 'default'");
    failures += report(load_recipe(lib_path, "").to_json() == Recipe{}.to_json(), "empty name loads 'default'");

    nlohmann::json lib;
    std::ifstream(lib_path) >> lib;
    int bad = 0;
    for (const auto &el : lib["roles"].items()) bad += throws_runtime([&] { load_recipe(lib_path, el.key()); });
    failures += report(bad == 0, "every role loads (no unknown keys)", lib["roles"].size());
    failures += report(throws_runtime([&] { load_recipe(lib_path, "best"); }), "a portfolio is refused");
    failures += report(throws_runtime([&] { load_recipe(lib_path, "no_such_recipe"); }), "an unknown name is refused");

    // a save_recipe-style file: an old problem key is skipped, an unknown knob is an error
    const auto tmp = std::filesystem::temp_directory_path() / "gropt_test_recipe.json";
    nlohmann::json entry = Recipe{}.to_json();
    entry["diff"]["concomitant_project"] = false;
    entry["diff"]["w_pns"] = 3.0;
    std::ofstream(tmp) << nlohmann::json{{"old", entry}}.dump();
    failures += report(load_recipe(tmp.string(), "").w_pns == 3.0, "old problem key skipped, knobs applied");
    entry["diff"].erase("safe_lifted");
    std::ofstream(tmp) << nlohmann::json{{"pre_lift", entry}}.dump();
    failures += report(!load_recipe(tmp.string(), "").safe_lifted, "a recipe without safe_lifted runs unlifted");
    entry["solver"]["no_such_knob"] = 1;
    std::ofstream(tmp) << nlohmann::json{{"bad", entry}}.dump();
    failures += report(throws_runtime([&] { load_recipe(tmp.string(), "bad"); }), "an unknown knob is an error");
    std::filesystem::remove(tmp);
    return failures;
}
