#ifndef WORKSPACE_SDMM_H
#define WORKSPACE_SDMM_H

#include "workspace_solver.hpp"

namespace Gropt {

class Operator;  // Forward declaration

struct WorkspaceSDMM : WorkspaceSolver {
    // SDMM-specific reweighting history
    Eigen::VectorXd yhat1, dyhat, dy, dhhat, dghat;
    Eigen::VectorXd yhat00, y00, s00, z00;

    void init(int Ax_size);
    void reinit(int Ax_size);
    void prep(Operator &op, const Eigen::VectorXd &X);
    void reweight(double rw_eps, double e_corr, double rw_scalelim);
};

}  // namespace Gropt

#endif
