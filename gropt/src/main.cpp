#include "spdlog/spdlog.h"

#include <exception>
#include <string>

void demo_diffusion(const std::string &recipe_path, const std::string &recipe_name);  // demo_diffusion.cpp

// Usage: gropt [recipe.json [recipe_name]]
//   no args           -> the built-in default recipe
//   recipe.json       -> the first recipe in that library file
//   recipe.json NAME  -> that named recipe
int main(int argc, char **argv){
    spdlog::set_level(spdlog::level::debug);
    const std::string recipe_path = (argc > 1) ? argv[1] : "";
    const std::string recipe_name = (argc > 2) ? argv[2] : "";
    try {
        demo_diffusion(recipe_path, recipe_name);
    } catch (const std::exception &e) {   // unreadable recipe file / unknown recipe name / bad mode
        spdlog::error("{}", e.what());
        return 1;
    }
    return 0;
}
