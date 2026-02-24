#ifndef OSQP_WORKSPACE_H
#define OSQP_WORKSPACE_H

#include "solver_workspace.hpp"

namespace Gropt {

// OSQPWorkspace uses the same per-operator state as the base SolverWorkspace.
// It can be extended here if OSQP-specific state is needed in the future.
using OSQPWorkspace = SolverWorkspace;

}  // namespace Gropt

#endif
