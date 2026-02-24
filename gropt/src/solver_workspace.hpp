#ifndef SOLVER_WORKSPACE_H
#define SOLVER_WORKSPACE_H

#include "Eigen/Dense"

namespace Gropt {

class Operator;  // Forward declaration

struct SolverWorkspace {
    // Per-operator iteration state (shared by all solvers)
    double weight = 1.0;
    double gamma = 1.5;

    bool do_rw = true;
    bool do_gamma = true;
    bool do_weight = true;
    bool do_scalelim = true;

    Eigen::VectorXd y0, y1;
    Eigen::VectorXd z0, z1;
    Eigen::VectorXd s0, s1;

    void init(int Ax_size);
    void reinit(int Ax_size);
    void prep(Operator &op, const Eigen::VectorXd &X);
};

}  // namespace Gropt

#endif
