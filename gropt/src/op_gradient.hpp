#ifndef OP_GRADIENT_H
#define OP_GRADIENT_H

/**
 * Constraint on gradient amplitude.  Supports the 'rot_variant' variable
 * to decide if gmax operates per axis or on the gradient magnitude
 *
 * Checks 'set_vals' and forces those values if they are not NaN
 */

#include <iostream>
#include <string>
#include <math.h>
#include "Eigen/Dense"

#include "op_main.hpp"

namespace Gropt {

class Op_Gradient : public Operator
{
    protected:
        double gmax;
        Eigen::VectorXd gmax_vec;

    public:
        Op_Gradient(const ProblemData &_pdata, double _gmax, bool _rot_variant, double _weight_mod);
        Op_Gradient(const ProblemData &_pdata, const Eigen::VectorXd &_gmax_vec, bool _rot_variant, double _weight_mod);
        virtual void init();

        virtual void forward(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void transpose(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void prox(Eigen::VectorXd &X);
        virtual void check(Eigen::VectorXd &X);
};

}  // close "namespace Gropt"

#endif
