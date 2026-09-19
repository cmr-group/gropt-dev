#ifndef OP_GRADIENT_H
#define OP_GRADIENT_H

/**
 * Gradient amplitude constraint |g| <= gmax [T/m]: per axis (rot_variant) or on the cross-axis magnitude.
 * The prox also pins samples with a non-NaN set_vals to that value.
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
