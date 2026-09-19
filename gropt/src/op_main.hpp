#ifndef OP_MAIN_H
#define OP_MAIN_H

/**
 * Base class for every constraint, regularization and objective term in GrOpt.
 * Subclasses implement forward, transpose and prox; the base class applies the spec_norm and
 * equilibration scaling and adds the ADMM and objective terms to the x-subproblem.
 */

#include "Eigen/Dense"
#include <iostream>
#include <math.h>
#include <string>
#include <vector>

#include "problem_data.hpp"

namespace Gropt {

struct WorkspaceSolver;

class Operator
{
  public:
    std::string name;

    // Warm-start key "<name>#<occurrence>", assigned in GroptParams::prepare()
    std::string unique_name;

    const ProblemData *pdata;

    int N;
    int Naxis;
    int Ntot;
    double dt;
    int Ax_size;

    bool rot_variant = true;
    bool do_init_weights = true;

    double target;
    double tol0;
    double tol;
    double cushion = 1e-2; // prox uses tol = (1 - cushion) * tol0, slightly inside the check() limit

    double spec_norm2;
    double spec_norm;

    double obj_weight = 1.0;
    // Objective feasibility gate in [0,1] (1 = off); set by the solver when obj_gate_enable
    double obj_gate = 1.0;
    double weight_mod = 1.0;
    double admm_weight = 1.0;
    // Objective ops only: add the objective as a gradient frozen at x0 in the RHS (DCA, via add_obj_rhs)
    // instead of as curvature in the LHS (add_obj); used for concave maximization terms
    bool linearize_obj = false;

    bool use_projection = false;

    Eigen::VectorXd x_temp;
    Eigen::VectorXd Ax_temp;

    bool do_equil = false;
    Eigen::VectorXd eq_rows;
    Eigen::VectorXd eq_cols;

    double feas_check;
    double r_feas;
    Eigen::VectorXd feas_temp;

    std::vector<int> hist_feas;
    std::vector<double> hist_r_feas;

    // ----------------------------------------
    Operator() = default;
    Operator(const ProblemData &_pdata);
    virtual ~Operator();

    virtual void init();

    virtual void forward(Eigen::VectorXd &X, Eigen::VectorXd &out);
    virtual void forward_op(Eigen::VectorXd &X, Eigen::VectorXd &out);
    virtual void transpose(Eigen::VectorXd &X, Eigen::VectorXd &out);
    virtual void transpose_op(Eigen::VectorXd &X, Eigen::VectorXd &out, bool apply_fixer);
    virtual void transpose_op(Eigen::VectorXd &X, Eigen::VectorXd &out);
    virtual void check(Eigen::VectorXd &X);
    virtual void prox(Eigen::VectorXd &X);
    virtual void get_feas(Eigen::VectorXd &s);
    virtual void add_Atb(Eigen::VectorXd &b, const WorkspaceSolver &ws);
    void add_AtAx(Eigen::VectorXd &x, Eigen::VectorXd &out, const WorkspaceSolver &ws);

    // Estimate this operator's spectral norm ||A|| by power-iterating its own AᵀA
    double estimate_self_spec_norm(int n_iters = 30);
    virtual void add_obj(Eigen::VectorXd &x, Eigen::VectorXd &out);
    virtual void add_obj_rhs(Eigen::VectorXd &x0, Eigen::VectorXd &out, bool normalize = false);

    // Freeze/unfreeze a nonlinear operator's linearization at X for the inner CG solve (no-op if linear)
    virtual void freeze_linearization(Eigen::VectorXd &X) {}
    virtual void unfreeze_linearization() {}

    // ||true(x_new) - frozen_linear(x_new)|| / ||true(x_new)||, 0 for linear ops; call while frozen.
    // Used by LinearizationErrorMonitor to reject steps that leave a nonlinear op's valid region.
    virtual double linearization_error(const Eigen::VectorXd &x_new) { return 0.0; }

    // Worst-sample true (nonlinear) constraint overage at x_new, 0 when feasible; only Op_SAFE implements it.
    // Used by FeasibilityMonitor and the objective gate.
    virtual double constraint_violation(const Eigen::VectorXd &x_new) { return 0.0; }
    void print_details();

    // EqualityProjection: append this op's equality row(s) and target(s), linearized at x0 if nonlinear;
    // nothing unless use_projection.
    virtual void append_eq_rows(std::vector<Eigen::VectorXd> &rows, std::vector<double> &targets,
                                const Eigen::VectorXd &x0) const {}
    // True if the appended rows depend on x0 (projector is rebuilt each outer iteration).
    virtual bool eq_rows_vary() const { return false; }

    // Partition of the Ax-space output (sums to Ax_size); warm start resizes the dual y block by block.
    virtual std::vector<int> Ax_block_lengths() const;
};

} // namespace Gropt

#endif
