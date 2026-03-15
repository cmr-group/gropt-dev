#ifndef OP_EDDY_H
#define OP_EDDY_H

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

class Op_Eddy : public Operator {

  protected:
    Eigen::MatrixXd A;
    Eigen::VectorXd lambdas;
    int Nlam;
    int Nrows;

  public:
    Op_Eddy(const ProblemData &_pdata, const Eigen::VectorXd &_lam, double _tol, double _weight_mod);

    virtual void init();

    virtual void forward(Eigen::VectorXd &X, Eigen::VectorXd &out);
    virtual void transpose(Eigen::VectorXd &X, Eigen::VectorXd &out);
    virtual void prox(Eigen::VectorXd &X);
    virtual void check(Eigen::VectorXd &X);
};

} // namespace Gropt

#endif
