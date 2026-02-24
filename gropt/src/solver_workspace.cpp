#include "solver_workspace.hpp"
#include "op_main.hpp"

namespace Gropt {

void SolverWorkspace::init(int Ax_size)
{
    y0.setZero(Ax_size);
    y1.setZero(Ax_size);
    z0.setZero(Ax_size);
    z1.setZero(Ax_size);
    s0.setZero(Ax_size);
    s1.setZero(Ax_size);
}

void SolverWorkspace::reinit(int Ax_size)
{
    if (y0.size() != Ax_size) { y0.setZero(Ax_size); }
    if (y1.size() != Ax_size) { y1.setZero(Ax_size); }
    if (s0.size() != Ax_size) { s0.setZero(Ax_size); }
    if (s1.size() != Ax_size) { s1.setZero(Ax_size); }
    z0.setZero(Ax_size);
    z1.setZero(Ax_size);
}

void SolverWorkspace::prep(Operator &op, const Eigen::VectorXd &X)
{
    Eigen::VectorXd X_copy = X;  // forward_op takes non-const ref
    op.forward_op(X_copy, z0);
    z1 = z0;
}

}  // namespace Gropt
