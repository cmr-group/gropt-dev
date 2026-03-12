#ifndef WORKSPACE_OSQP_H
#define WORKSPACE_OSQP_H

#include "workspace_solver.hpp"

namespace Gropt {

struct WorkspaceOSQP : WorkspaceSolver {
    Eigen::VectorXd y00, z00;
    double weight_scale = 1.0;

    void init(int Ax_size);
    void reinit(int Ax_size);
    void prep(Operator &op, const Eigen::VectorXd &X);
    void init_reweight();
    void reweight(int iiter, double scale_up, double scale_down, double eps);
};

} // namespace Gropt

#endif
