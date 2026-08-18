#ifndef OP_MAIN_H
#define OP_MAIN_H

/**
 * This is our main parent class for every constraint and regularization term in GrOpt
 * The main functions that we use here are the reweighting schemes, where most operators will
 * then just need to implement the forward, transpose, and proximal mapping operators.
 */

#include "Eigen/Dense"
#include <iostream>
#include <math.h>
#include <string>
#include <vector>

#include "problem_data.hpp"

namespace Gropt {

struct WorkspaceSolver; // Forward declaration

class Operator // This is the main parent class for every operator in GrOpt
{
  public:
    std::string name;

    // Stable per-operator identifier, assigned in GroptParams::prepare() as "<name>#<occurrence>"
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
    double cushion = 1e-2; // This is the cushion factor, which is used to reduce the tolerance by a factor of 1-cushion

    double spec_norm2;
    double spec_norm;

    double obj_weight = 1.0;
    // Feasibility gate in [0,1] applied to the linearized objective pull (add_obj_rhs). The solver sets
    // it each outer iteration from the total constraint violation (~0 while infeasible, ->1 when
    // feasible) so the objective can't cold-launch past the constraints before they engage; 1.0 = off.
    double obj_gate = 1.0;
    double weight_mod = 1.0;
    double admm_weight = 1.0;
    // Objective ops only: if true, this objective enters the x-subproblem as a linearized gradient
    // in the RHS 
    bool linearize_obj = false;
    // Objective ops only: true if this is the optimization TARGET used for best-feasible scoring
    // (e.g. b-value max). False for objective-path terms that merely drive a constraint to zero
    bool is_score_obj = true;

    bool use_projection = false;

    bool fix_gamma = false;
    double gamma_fix = 1.0;

    Eigen::VectorXd x_temp;
    Eigen::VectorXd x_temp_obj;
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
    // Refresh per-outer-iteration objective state at the current iterate x0 (default no-op). Used by
    // augmented-Lagrangian objective terms to freeze their linearization and update the scalar dual.
    virtual void update_obj_state(const Eigen::VectorXd &x0) {}

    // Freeze/unfreeze the operator's linearization for the duration of the inner CG solve (default
    // no-op; linear ops need nothing).
    virtual void freeze_linearization(Eigen::VectorXd &X) {}
    virtual void unfreeze_linearization() {}

    // Model-fidelity of the frozen linearization at a candidate iterate x_new: relative error between the
    // TRUE (nonlinear) forward and the FROZEN-LINEAR forward, ||true - linear|| / ||true||. 0 for linear
    // operators (their linearization is exact). Used by the trust-region StepMonitor to reject a step that
    // outran a nonlinear operator's (e.g. SAFE's) valid linearization region. Assumes the operator is
    // still frozen at the outer iterate (called between freeze_linearization and unfreeze_linearization).
    virtual double linearization_error(const Eigen::VectorXd &x_new) { return 0.0; }

    // Total TRUE (nonlinear) constraint violation at x_new: sum over the operator's rows of the amount
    // each exceeds its limit (0 when feasible, and 0 for objectives / non-constraint ops). Used by the
    // merit-based StepMonitor to penalize stepping past the constraint. Independent of the frozen
    // linearization (evaluates the real constraint).
    virtual double constraint_violation(const Eigen::VectorXd &x_new) { return 0.0; }
    void print_details();

    // EqualityProjection participation. An operator enforcing a linear equality (or a nonlinear one
    // relinearized at the current iterate x0) appends its coefficient row(s) and target(s) here.
    // Operators not using projection append nothing.
    virtual void append_eq_rows(std::vector<Eigen::VectorXd> &rows, std::vector<double> &targets,
                                const Eigen::VectorXd &x0) const {}
    // True if the appended rows depend on x0 (relinearized each outer iteration -> rebuild projector).
    virtual bool eq_rows_vary() const { return false; }

    // Flat partition of this operator's Ax-space output (sum == Ax_size), used to resize the dual
    // y across a grid change one block at a time so stacked blocks never bleed together (see
    // warmstart.hpp).
    virtual std::vector<int> Ax_block_lengths() const;
};

} // namespace Gropt

#endif
