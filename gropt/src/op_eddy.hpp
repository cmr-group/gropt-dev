#ifndef OP_EDDY_H
#define OP_EDDY_H

/**
 * Residual eddy-current constraint: for each axis and time constant lambda [s], the eddy current left at
 * the end of the waveform (exponential-kernel functional of g) must satisfy |e| <= tol (target 0).
 * With use_projection it is instead driven to exactly 0 by the equality projection.
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
    virtual void append_eq_rows(std::vector<Eigen::VectorXd> &rows, std::vector<double> &targets,
                                const Eigen::VectorXd &x0) const;
};

} // namespace Gropt

#endif
