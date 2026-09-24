#include <stdexcept>

#include "problem_data.hpp"

namespace Gropt {

int ProblemData::add_aux(const std::string &name, int size) {
    for (const auto &a : aux) {
        if (a.name == name) {
            if (a.size != size) {
                throw std::invalid_argument("ProblemData::add_aux: block '" + name + "' already exists with "
                                            "size " + std::to_string(a.size) + ", cannot redeclare it as " +
                                            std::to_string(size));
            }
            return a.offset;
        }
    }
    if (size <= 0) throw std::invalid_argument("ProblemData::add_aux: size must be > 0");
    const int offset = n_total();
    aux.push_back({name, offset, size});
    return offset;
}

} // namespace Gropt
