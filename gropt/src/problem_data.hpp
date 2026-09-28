#ifndef PROBLEM_DATA_H
#define PROBLEM_DATA_H

#include <stdexcept>
#include <string>
#include <vector>

#include "Eigen/Dense"

namespace Gropt {

// Auxiliary primal block after the waveform (e.g. a lifted SAFE's u >= |slew|); registered by name so
// operators needing the same quantity share it.
struct AuxBlock {
    std::string name;
    int offset = 0; // index of the block's first entry in the full primal vector
    int size = 0;
};

struct ProblemData {
    int N = 10;
    int Naxis = 1;
    double dt = 10e-6;

    Eigen::VectorXd X0;       // Initial guess
    Eigen::VectorXd inv_vec;  // Inversion vector (diffusion encoding sign flips)
    Eigen::VectorXd set_vals; // Fixed values (NaN = free); waveform only
    Eigen::VectorXd fixer;    // Binary mask over the FULL primal: 0 = fixed, 1 = free

    // Auxiliary blocks, in registration order. Rebuilt every prepare().
    std::vector<AuxBlock> aux;

    int n_wave() const { return N * Naxis; }
    int n_total() const {
        int n = n_wave();
        for (const auto &a : aux) n += a.size;
        return n;
    }

    // Offset of a registered block, or -1 if it was never declared.
    int aux_offset(const std::string &name) const {
        for (const auto &a : aux) {
            if (a.name == name) return a.offset;
        }
        return -1;
    }

    // Register a block, or return the existing one of the same name (redeclaring another size throws).
    int add_aux(const std::string &name, int size) {
        for (const auto &a : aux) {
            if (a.name != name) continue;
            if (a.size != size) {
                throw std::invalid_argument("ProblemData::add_aux: block '" + name + "' already exists with size " +
                                            std::to_string(a.size) + ", cannot redeclare it as " + std::to_string(size));
            }
            return a.offset;
        }
        if (size <= 0) throw std::invalid_argument("ProblemData::add_aux: size must be > 0");
        aux.push_back({name, n_total(), size});
        return aux.back().offset;
    }

    void clear_aux() { aux.clear(); }
};

}  // namespace Gropt

#endif
