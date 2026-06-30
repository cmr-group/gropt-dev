#ifndef OP_CONCOMITANT_H
#define OP_CONCOMITANT_H

/**
 * Constraint on gradient amplitude.  Supports the 'rot_variant' variable
 * to decide if gmax operates per axis or on the gradient magnitude
 *
 * Checks 'set_vals' and forces those values if they are not NaN
 */

#include "Eigen/Dense"
#include <iostream>
#include <math.h>
#include <string>

#include "op_main.hpp"

namespace Gropt {

class Op_Concomitant : public Operator {

  protected:
    int start_idx = 0;
    double con_tol0 = 0.1;     // true fractional tolerance on the pre/post energy ratio
    double con_target = 1.0;   // target pre/post energy ratio pos/neg (1 = balanced)

    // Augmented-Lagrangian state, used only in the objective (all_obj) mode. The nonconvex equality
    // c(x)=pos-neg=0 is driven via lambda*c + (rho/2)c^2, linearized at the iterate each outer iter.
    double auglag_rho = 1.0;     // penalty parameter (= weight_mod)
    double auglag_lambda = 0.0;  // scalar dual (method of multipliers)
    double auglag_c = 0.0;       // frozen c(x0) = pos - neg
    double auglag_gx0 = 0.0;     // frozen g . x0
    Eigen::VectorXd auglag_g;    // frozen masked gradient grad c(x0)

  public:
    Op_Concomitant(const ProblemData &_pdata, int _start_idx, bool _rot_variant, double _weight_mod,
                   double _tol0 = 0.1, double _target = 1.0);
    virtual void init();

    virtual void forward(Eigen::VectorXd &X, Eigen::VectorXd &out);
    virtual void transpose(Eigen::VectorXd &X, Eigen::VectorXd &out);
    virtual void prox(Eigen::VectorXd &X);
    virtual void check(Eigen::VectorXd &X);

    // Nonconvex energy-balance equality c(x)=pos-neg=0, linearized at x0 (SQP): one row
    // a = grad c(x0), target c0. eq_rows_vary is true (relinearized each outer iteration).
    virtual void append_eq_rows(std::vector<Eigen::VectorXd> &rows, std::vector<double> &targets,
                                const Eigen::VectorXd &x0) const;
    virtual bool eq_rows_vary() const { return true; }

    // Objective (all_obj) mode: augmented-Lagrangian handling of c(x)=0. update_obj_state freezes
    // the linearization (g, c0) at x0 and advances the scalar dual; add_obj/add_obj_rhs inject the
    // rank-1 Gauss-Newton curvature (LHS) and the dual/penalty pull (RHS).
    virtual void update_obj_state(const Eigen::VectorXd &x0);
    virtual void add_obj(Eigen::VectorXd &X, Eigen::VectorXd &out);
    virtual void add_obj_rhs(Eigen::VectorXd &x0, Eigen::VectorXd &out, bool normalize);
};

} // namespace Gropt

#endif
