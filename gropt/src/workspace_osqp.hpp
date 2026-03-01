#ifndef WORKSPACE_OSQP_H
#define WORKSPACE_OSQP_H

#include "workspace_solver.hpp"

namespace Gropt {

// WorkspaceOSQP uses the same per-operator state as the base WorkspaceSolver.
// It can be extended here if OSQP-specific state is needed in the future.
using WorkspaceOSQP = WorkspaceSolver;

}  // namespace Gropt

#endif
