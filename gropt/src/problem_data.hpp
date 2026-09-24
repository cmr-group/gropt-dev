#ifndef PROBLEM_DATA_H
#define PROBLEM_DATA_H

#include <string>
#include <vector>

#include "Eigen/Dense"

namespace Gropt {

// An auxiliary primal block appended after the waveform. Lifting a nonsmooth term (e.g. SAFE's
// LP2(|slew|), via u >= |slew|) needs a variable the optimizer chooses, which cannot live inside an
// operator: operators are maps FROM the primal. Blocks are registered BY NAME so that several operators
// needing the same quantity share one -- the PNS and cardiac Op_SAFE both read one "abs_slew" block --
// while a future constraint needing something else registers its own.
struct AuxBlock {
    std::string name;
    int offset = 0; // index of the block's first entry in the full primal vector
    int size = 0;
};

struct ProblemData {
    int N = 10;
    int Naxis = 1;
    double dt = 10e-6;

    Eigen::VectorXd X0;       // Initial guess (waveform only; aux blocks are seeded by their operator)
    Eigen::VectorXd inv_vec;  // Inversion vector (diffusion encoding sign flips)
    Eigen::VectorXd set_vals; // Fixed values (NaN = free); waveform only
    Eigen::VectorXd fixer;    // Binary mask over the FULL primal: 0 = fixed, 1 = free

    // Auxiliary blocks, in registration order. Rebuilt every prepare().
    std::vector<AuxBlock> aux;

    int n_wave() const { return N * Naxis; }

    int n_aux() const {
        int n = 0;
        for (const auto &a : aux) n += a.size;
        return n;
    }

    int n_total() const { return n_wave() + n_aux(); }

    // Offset of a registered block, or -1 if it was never declared.
    int aux_offset(const std::string &name) const {
        for (const auto &a : aux) {
            if (a.name == name) return a.offset;
        }
        return -1;
    }

    int aux_size(const std::string &name) const {
        for (const auto &a : aux) {
            if (a.name == name) return a.size;
        }
        return 0;
    }

    // Register a block, or return the existing one of the same name. Declaring the same name with a
    // different size is a programming error: the two operators disagree about what they are sharing.
    int add_aux(const std::string &name, int size);

    void clear_aux() { aux.clear(); }
};

} // namespace Gropt

#endif
