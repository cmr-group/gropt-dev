#ifndef OP_CONCOMITANT_H
#define OP_CONCOMITANT_H

/**
 * Concomitant constraint: balances the gradient energy sum g^2 dt before (pos) and after (neg) the 180,
 * |pos/neg - target| <= tol0. Samples before start_idx are ignored.
 * exact_quad (default) uses the exact piecewise-linear integral sum (dt/3)(g_i^2 + g_i g_{i+1} + g_{i+1}^2),
 * so the ratio is raster independent; the legacy sum is low by (dt^2/6) * integral(slew^2) per side.
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
    double con_tol0 = 0.1;     // tolerance band on pos/neg
    double con_target = 1.0;   // target pre/post energy ratio pos/neg (1 = balanced)

    // Pre/post energies of X, and d(pos)/dX and d(neg)/dX (either gradient may be null if not wanted).
    void energies(const Eigen::VectorXd &X, double &pos, double &neg, Eigen::VectorXd *d_pos,
                  Eigen::VectorXd *d_neg) const;

  public:
    Op_Concomitant(const ProblemData &_pdata, int _start_idx, bool _rot_variant, double _weight_mod,
                   double _tol0 = 0.1, double _target = 1.0);
    virtual void init();

    // false = legacy rectangle sum g^2 dt, kept for comparison
    bool exact_quad = true;

    virtual void forward(Eigen::VectorXd &X, Eigen::VectorXd &out);
    virtual void transpose(Eigen::VectorXd &X, Eigen::VectorXd &out);
    virtual void prox(Eigen::VectorXd &X);
    virtual void check(Eigen::VectorXd &X);

    // Nonconvex equality c(x) = pos - target*neg = 0 linearized at x0: row grad c(x0), rhs c(x0)
    virtual void append_eq_rows(std::vector<Eigen::VectorXd> &rows, std::vector<double> &targets,
                                const Eigen::VectorXd &x0) const;
    virtual bool eq_rows_vary() const { return true; }
};

} // namespace Gropt

#endif
